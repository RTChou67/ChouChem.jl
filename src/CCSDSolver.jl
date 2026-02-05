export spin_orbital_ccsd, CCSDResults

struct CCSDResults
	EtotHF::Float64
	EtotCCSD::Float64
	Ecorr::Float64
	T1::Matrix{Float64}
	T2::Array{Float64, 4}
end

function spin_orbital_ccsd(SCF_Res::UHFResults; MaxIter=50, Threshold=1e-8)
	println("\n--- Starting CCSD Calculation (Spin-Orbital) ---")
	
	ONum = length(SCF_Res.BasisSet)
	ENum = SCF_Res.ENumAlpha + SCF_Res.ENumBeta
	SONum = 2 * ONum
	
	C_a, C_b = SCF_Res.CAlpha, SCF_Res.CBeta
	H_a = C_a' * SCF_Res.Hcore * C_a
	H_b = C_b' * SCF_Res.Hcore * C_b
	
	eps_a = SCF_Res.EAlpha
	eps_b = SCF_Res.EBeta
	
	N_occ_a = SCF_Res.ENumAlpha
	N_occ_b = SCF_Res.ENumBeta
	N_vir_a = ONum - N_occ_a
	N_vir_b = ONum - N_occ_b
	
	No = N_occ_a + N_occ_b
	Nv = N_vir_a + N_vir_b
	
	C_spin = zeros(Float64, ONum, SONum) 
	
	println("Transforming integrals...")
	@time begin
		ERI_aa = TransERI_blas(SCF_Res.ERI, C_a, C_a, C_a, C_a, ONum)
		ERI_bb = TransERI_blas(SCF_Res.ERI, C_b, C_b, C_b, C_b, ONum)
		ERI_ab = TransERI_blas(SCF_Res.ERI, C_a, C_a, C_b, C_b, ONum)
		
		V_so_raw = build_spin_eri(ERI_aa, ERI_bb, ERI_ab, ONum)
	end

	V_so_phys = permutedims(V_so_raw, (1, 3, 2, 4))
	
	perm = Vector{Int}()
	append!(perm, 1:N_occ_a)
	append!(perm, (ONum+1):(ONum+N_occ_b))
	append!(perm, (N_occ_a+1):ONum)
	append!(perm, (ONum+N_occ_b+1):(2*ONum))
	
V = V_so_phys[perm, perm, perm, perm]
	
f = zeros(Float64, SONum)
	f[1:N_occ_a] = eps_a[1:N_occ_a]
	f[N_occ_a+1:No] = eps_b[1:N_occ_b]
	f[No+1:No+N_vir_a] = eps_a[N_occ_a+1:end]
	f[No+N_vir_a+1:end] = eps_b[N_occ_b+1:end]
	
V_oooo = V[1:No, 1:No, 1:No, 1:No]
	V_vvoo = V[No+1:end, No+1:end, 1:No, 1:No]
	V_voov = V[No+1:end, 1:No, 1:No, No+1:end] 
    
    get_v(p, q, r, s) = V[p, q, r, s] - V[p, q, s, r]
    
    T1 = zeros(Float64, No, Nv)
    T2 = zeros(Float64, No, No, Nv, Nv)
    
    D1 = [f[i] - f[No+a] for i in 1:No, a in 1:Nv]
    D2 = [f[i] + f[j] - f[No+a] - f[No+b] for i in 1:No, j in 1:No, a in 1:Nv, b in 1:Nv]
    
    for i in 1:No, j in 1:No, a in 1:Nv, b in 1:Nv
        T2[i,j,a,b] = get_v(i, j, No+a, No+b) / D2[i,j,a,b]
    end
    
    DIIS = DIISManager{Tuple{Matrix{Float64}, Array{Float64,4}}}(8)
    
    E_corr_old = 0.0
    println("\nIter     E_corr        Delta_E       RMS(Res)")
    println("----------------------------------------------")
    
    for iter in 1:MaxIter
        Tau = zeros(Float64, No, No, Nv, Nv) 
        Tau_tilde = zeros(Float64, No, No, Nv, Nv) 
        
        @inbounds for b in 1:Nv, a in 1:Nv, j in 1:No, i in 1:No
            p = T1[i,a] * T1[j,b]
            t2 = T2[i,j,a,b]
            Tau[i,j,a,b] = t2 + p - T1[i,b] * T1[j,a] 
            Tau_tilde[i,j,a,b] = t2 + 0.5 * (p - T1[i,b] * T1[j,a])
        end
        
        F_ae = zeros(Float64, Nv, Nv)
        for e in 1:Nv, a in 1:Nv
            val = (a==e ? f[No+a] : 0.0)
            sum_mf = 0.0
            for m in 1:No, f in 1:Nv
                sum_mf += T1[m,f] * get_v(m, No+a, No+f, No+e)
            end
            val += sum_mf
            sum_mef = 0.0
            for m in 1:No, n in 1:No, f in 1:Nv
                 sum_mef += Tau_tilde[m,n,a,f] * get_v(m, n, No+e, No+f)
            end
            val -= 0.5 * sum_mef
            F_ae[a,e] = val
        end
        
        F_mi = zeros(Float64, No, No)
        for i in 1:No, m in 1:No
            val = (m==i ? f[m] : 0.0)
            sum_ne = 0.0
            for n in 1:No, e in 1:Nv
                sum_ne += T1[n,e] * get_v(m, n, i, No+e)
            end
            val += sum_ne
            sum_nef = 0.0
            for n in 1:No, e in 1:Nv, f in 1:Nv
                sum_nef += Tau_tilde[i,n,e,f] * get_v(m, n, No+e, No+f)
            end
            val += 0.5 * sum_nef
            F_mi[m,i] = val
        end
        
        F_me = zeros(Float64, No, Nv)
        for e in 1:Nv, m in 1:No
            val = 0.0 
            for n in 1:No, f in 1:Nv
                val += T1[n,f] * get_v(m, n, No+e, No+f)
            end
            F_me[m,e] = val
        end
        
        W_mnij = zeros(Float64, No, No, No, No)
        for j in 1:No, i in 1:No, n in 1:No, m in 1:No
            val = get_v(m, n, i, j)
            val += sum(T1[j,e] * get_v(m, n, i, No+e) for e in 1:Nv)
            val -= sum(T1[i,e] * get_v(m, n, j, No+e) for e in 1:Nv)
            val += 0.25 * sum(Tau[i,j,e,f] * get_v(m, n, No+e, No+f) for e in 1:Nv, f in 1:Nv)
            W_mnij[m,n,i,j] = val
        end
        
        W_abef = zeros(Float64, Nv, Nv, Nv, Nv)
        for f in 1:Nv, e in 1:Nv, b in 1:Nv, a in 1:Nv
            val = get_v(No+a, No+b, No+e, No+f)
            val -= sum(T1[m,b] * get_v(m, No+a, No+e, No+f) for m in 1:No)
            val += sum(T1[m,a] * get_v(m, No+b, No+e, No+f) for m in 1:No)
            val += 0.25 * sum(Tau[m,n,a,b] * get_v(m, n, No+e, No+f) for m in 1:No, n in 1:No)
            W_abef[a,b,e,f] = val
        end
        
        W_mbej = zeros(Float64, No, Nv, Nv, No)
        for j in 1:No, e in 1:Nv, b in 1:Nv, m in 1:No
            val = get_v(m, No+b, No+e, j)
            val += sum(T1[j,f] * get_v(m, No+b, No+e, No+f) for f in 1:Nv)
            val -= sum(T1[n,b] * get_v(m, n, No+e, j) for n in 1:No)
            val -= sum( (0.5 * T2[j,n,f,b] + T1[j,f]*T1[n,b]) * get_v(m, n, No+e, No+f) for n in 1:No, f in 1:Nv)
            W_mbej[m,b,e,j] = val
        end
        
        R1 = zeros(Float64, No, Nv)
        for a in 1:Nv, i in 1:No
            val = 0.0 
            val += sum(T1[i,e] * F_ae[a,e] for e in 1:Nv)
            val -= sum(T1[m,a] * F_mi[m,i] for m in 1:No)
            val += sum(T2[i,m,a,e] * F_me[m,e] for m in 1:No, e in 1:Nv)
            val -= sum(T1[n,f] * get_v(n, No+a, i, No+f) for n in 1:No, f in 1:Nv)
            val -= 0.5 * sum(T2[i,m,e,f] * get_v(m, No+a, No+e, No+f) for m in 1:No, e in 1:Nv, f in 1:Nv)
            val -= 0.5 * sum(T2[m,n,a,e] * get_v(n, m, i, No+e) for m in 1:No, n in 1:No, e in 1:Nv)
            R1[i,a] = val
        end
        
        R2 = zeros(Float64, No, No, Nv, Nv)
        for b in 1:Nv, a in 1:Nv, j in 1:No, i in 1:No
            val = get_v(i, j, No+a, No+b)
            term1 = 0.0
            for e in 1:Nv
                term1 += T2[i,j,a,e] * (F_ae[b,e] - 0.5 * sum(T1[m,b] * F_me[m,e] for m in 1:No))
            end
            term1_perm = 0.0
            for e in 1:Nv
                term1_perm += T2[i,j,b,e] * (F_ae[a,e] - 0.5 * sum(T1[m,a] * F_me[m,e] for m in 1:No))
            end
            val += term1 - term1_perm
            
            term2 = 0.0
            for m in 1:No
                term2 += T2[i,m,a,b] * (F_mi[m,j] + 0.5 * sum(T1[j,e] * F_me[m,e] for e in 1:Nv))
            end
            term2_perm = 0.0
            for m in 1:No
                term2_perm += T2[j,m,a,b] * (F_mi[m,i] + 0.5 * sum(T1[i,e] * F_me[m,e] for e in 1:Nv))
            end
            val -= (term2 - term2_perm)
            
            val += 0.5 * sum(Tau[m,n,a,b] * W_mnij[m,n,i,j] for m in 1:No, n in 1:No)
            
            val += 0.5 * sum(Tau[i,j,e,f] * W_abef[a,b,e,f] for e in 1:Nv, f in 1:Nv)
            
            t_base = sum( (T2[i,m,a,e] * W_mbej[m,b,e,j] - T1[i,e]*T1[m,a]*get_v(m, No+b, No+e, j)) for m in 1:No, e in 1:Nv)
            t_p_ab = sum( (T2[i,m,b,e] * W_mbej[m,a,e,j] - T1[i,e]*T1[m,b]*get_v(m, No+a, No+e, j)) for m in 1:No, e in 1:Nv)
            t_p_ij = sum( (T2[j,m,a,e] * W_mbej[m,b,e,i] - T1[j,e]*T1[m,a]*get_v(m, No+b, No+e, i)) for m in 1:No, e in 1:Nv)
            t_p_all = sum( (T2[j,m,b,e] * W_mbej[m,a,e,i] - T1[j,e]*T1[m,b]*get_v(m, No+a, No+e, i)) for m in 1:No, e in 1:Nv)
            
            val += (t_base - t_p_ab - t_p_ij + t_p_all)
            
            val += sum(T1[i,e] * get_v(No+a, No+b, No+e, j) for e in 1:Nv)
            val -= sum(T1[j,e] * get_v(No+a, No+b, No+e, i) for e in 1:Nv)
            
            val -= sum(T1[m,a] * get_v(m, No+b, i, j) for m in 1:No)
            val += sum(T1[m,b] * get_v(m, No+a, i, j) for m in 1:No)
            
            R2[i,j,a,b] = val
        end
        
        E_corr = 0.0
        for i in 1:No, j in 1:No, a in 1:Nv, b in 1:Nv
            v_int = get_v(i, j, No+a, No+b)
            E_corr += 0.25 * v_int * T2[i,j,a,b] + 0.5 * v_int * T1[i,a] * T1[j,b]
        end
        
        delta_E = E_corr - E_corr_old
        E_corr_old = E_corr
        
        T1_new = T1 .+ R1 ./ D1
        T2_new = T2 .+ R2 ./ D2
        
        T1_final, T2_final = diis_update!(DIIS, (T1_new, T2_new), (R1, R2))
        
        rms_res = sqrt( (sum(R1.^2) + sum(R2.^2)) / (length(R1) + length(R2)) )
        
        @printf("%4d   %14.10f   %10.2e   %10.2e\n", iter, E_corr, delta_E, rms_res)
        
        T1 = T1_final
        T2 = T2_final
        
        if abs(delta_E) < Threshold && rms_res < Threshold
            println("\nCCSD converged.")
            return CCSDResults(SCF_Res.Etot, SCF_Res.Etot + E_corr, E_corr, T1, T2)
        end
    end
    
    println("\nCCSD failed to converge.")
    return nothing
end