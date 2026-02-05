using LinearAlgebra
using Printf
using ChemAlgebra: DIISManager, diis_update!

struct CCSDResults
	EtotHF::Float64
	EtotCCSD::Float64
	Ecorr::Float64
	T1::Matrix{Float64}
	T2::Array{Float64, 4}
end

function spin_orbital_ccsd(SCF_Res::UHFResults; MaxIter=50, Threshold=1e-8)
	println("\n--- Starting CCSD Calculation (Spin-Orbital) ---")
	
	# 1. Prepare Integrals and Fock Matrix in Spin-Orbitals
	ONum = length(SCF_Res.BasisSet)
	ENum = SCF_Res.ENumAlpha + SCF_Res.ENumBeta
	SONum = 2 * ONum
	
	C_a, C_b = SCF_Res.CAlpha, SCF_Res.CBeta
	H_a = C_a' * SCF_Res.Hcore * C_a
	H_b = C_b' * SCF_Res.Hcore * C_b
	
	# Fock Matrix (Spin-Orbital)
	# F_spin is block diagonal: [F_alpha 0; 0 F_beta]
	# However, we need the "true" Fock matrix elements including correlation potential if we weren't using canonical HF orbitals.
	# Since we use canonical HF orbitals, F is diagonal with orbital energies.
	# We construct the diagonal f_p^p vector directly from eigenvalues.
	eps_a = SCF_Res.EAlpha
	eps_b = SCF_Res.EBeta
	
	# Mapping: 1..ONum -> Alpha, ONum+1..2*ONum -> Beta
	# We need to sort them or handle them as is. Usually, we keep occ and vir blocks contiguous.
	# But standard CCSD codes often split into Occ and Vir ranges.
	# Let's define the Spin-Orbital ranges.
	# Since UHF results might sort orbitals by energy, we must ensure Occupied are first.
	# RHF/UHF usually return sorted eigenvalues. 
	
	# Occupied Indices: 1..N_alpha (Alpha), N_alpha+1..N_alpha+N_beta (Beta)?
	# No, the C matrices are usually N_basis x N_basis.
	# The first N_alpha columns of C_alpha are occupied.
	
	# Let's map spin orbitals to a single index p = 1..2M.
	# Order: Alpha(1..M), Beta(1..M).
	# Then we permute to: Occ(Alpha), Occ(Beta), Vir(Alpha), Vir(Beta).
	
	N_occ_a = SCF_Res.ENumAlpha
	N_occ_b = SCF_Res.ENumBeta
	N_vir_a = ONum - N_occ_a
	N_vir_b = ONum - N_occ_b
	
	No = N_occ_a + N_occ_b
	Nv = N_vir_a + N_vir_b
	
	# Construct MO coefficients in the order: [Occ_A, Occ_B, Vir_A, Vir_B]
	# This simplifies the loops: i,j,k,l in 1:No; a,b,c,d in 1:Nv.
	
	C_spin = zeros(Float64, ONum, SONum) # This structure is tricky. We need spatial AO to spin MO.
	# Actually, we transformation separately.
	
	println("Transforming integrals...")
	@time begin
		ERI_aa = TransERI_blas(SCF_Res.ERI, C_a, C_a, C_a, C_a, ONum)
		ERI_bb = TransERI_blas(SCF_Res.ERI, C_b, C_b, C_b, C_b, ONum)
		ERI_ab = TransERI_blas(SCF_Res.ERI, C_a, C_a, C_b, C_b, ONum)
		
		# V_so_raw indices: p, q, r, s in 1..2M (Alpha then Beta)
		V_so_raw = build_spin_eri(ERI_aa, ERI_bb, ERI_ab, ONum)
	end

	# Convert Chemist Notation (pq|rs) to Physics Notation <pr|qs>
	# Chemist: 1, 2, 3, 4 -> Electron 1 in (1,3), Electron 2 in (2,4) ??
	# Wait. (pq|rs) = \int p*(1) q(1) 1/r12 r*(2) s(2). 
	# Physics <pr|qs> = \int p*(1) r*(2) 1/r12 q(1) s(2) = \int p*(1) q(1) ... r*(2) s(2) = (pq|rs).
	# Mapping:
	# Chemist (p, q, r, s) -> (p q | r s).
	# Physics < p q | r s > = \int p*(1) q*(2) 1/r12 r(1) s(2)
	#                       = \int p*(1) r(1) 1/r12 q*(2) s(2)
	#                       = (p r | q s)
	# So V_phys[p, q, r, s] should be V_chem[p, r, q, s]
	
	V_so_phys = permutedims(V_so_raw, (1, 3, 2, 4))
	
	# Create permutation map to [Occ_A, Occ_B, Vir_A, Vir_B]
	# Current V_so_raw order: 1..M (Alpha), M+1..2M (Beta)
	# Alpha Occ: 1 .. N_occ_a
	# Alpha Vir: N_occ_a+1 .. M
	# Beta Occ: M+1 .. M+N_occ_b
	# Beta Vir: M+N_occ_b+1 .. 2M
	
	perm = Vector{Int}()
	# Occ
	append!(perm, 1:N_occ_a)
	append!(perm, (ONum+1):(ONum+N_occ_b))
	# Vir
	append!(perm, (N_occ_a+1):ONum)
	append!(perm, (ONum+N_occ_b+1):(2*ONum))
	
	# Permute V_so
	V = V_so_phys[perm, perm, perm, perm]
	
	# Construct Fock vector (diagonal) in the same order
	f = zeros(Float64, SONum)
	# Occ
	f[1:N_occ_a] = eps_a[1:N_occ_a]
	f[N_occ_a+1:No] = eps_b[1:N_occ_b]
	# Vir
	f[No+1:No+N_vir_a] = eps_a[N_occ_a+1:end]
	f[No+N_vir_a+1:end] = eps_b[N_occ_b+1:end]
	
	# Antisymmetrized Integrals: <ij||ab> = <ij|ab> - <ij|ba>
	# V in physicist notation is often used <pq|rs> -> (pr|qs).
	# build_spin_eri returns chemist notation <pq|rs> (1 2 | 1 2)
	# Physics notation <pq||rs> = <pq|rs> - <pq|sr>
	
	# Prepare tensor views
	# i,j,k,l,m,n -> 1:No
	# a,b,c,d,e,f -> 1:Nv (offset by No in global arrays if needed, but we split)
	
	V_oooo = V[1:No, 1:No, 1:No, 1:No]
	V_vvoo = V[No+1:end, No+1:end, 1:No, 1:No]
	V_voov = V[No+1:end, 1:No, 1:No, No+1:end] # <ai|jb> needed for T1
    # We need <mb||ej> -> <Occ, Vir || Vir, Occ>
    # In V array: V[m, b, e, j] -> V[1:No, No+1:end, No+1:end, 1:No]
    # Be careful with antisymmetrization!
    
    # Let's pre-calculate antisymmetrized integral slices needed
    # Function to get <pq||rs> = V[p,q,r,s] - V[p,q,s,r]
    get_v(p, q, r, s) = V[p, q, r, s] - V[p, q, s, r]
    
    # We need to construct full blocks of antisymmetrized integrals for efficiency
    # Or use loops. For "naive" implementation, loops with get_v is safest but slow.
    # We will build key blocks.
    
    # Initial Amplitudes (MP2 Guess)
    # D_ij^ab = f_i + f_j - f_a - f_b
    T1 = zeros(Float64, No, Nv)
    T2 = zeros(Float64, No, No, Nv, Nv)
    
    D1 = [f[i] - f[No+a] for i in 1:No, a in 1:Nv]
    D2 = [f[i] + f[j] - f[No+a] - f[No+b] for i in 1:No, j in 1:No, a in 1:Nv, b in 1:Nv]
    
    # MP2 T2 Guess
    for i in 1:No, j in 1:No, a in 1:Nv, b in 1:Nv
        T2[i,j,a,b] = get_v(i, j, No+a, No+b) / D2[i,j,a,b]
    end
    
    DIIS = DIISManager{Tuple{Matrix{Float64}, Array{Float64,4}}}(8)
    
    E_corr_old = 0.0
    println("\nIter     E_corr        Delta_E       RMS(Res)")
    println("----------------------------------------------")
    
    for iter in 1:MaxIter
        # Update Intermediates (Tau) 
        Tau = zeros(Float64, No, No, Nv, Nv) # t_ij^ab + t_i^a t_j^b
        Tau_tilde = zeros(Float64, No, No, Nv, Nv) # t_ij^ab + 0.5 t_i^a t_j^b
        
        # This loop is O(N^4), acceptable
        @inbounds for b in 1:Nv, a in 1:Nv, j in 1:No, i in 1:No
            p = T1[i,a] * T1[j,b]
            t2 = T2[i,j,a,b]
            Tau[i,j,a,b] = t2 + p - T1[i,b] * T1[j,a] # Antisymmetrized product
            Tau_tilde[i,j,a,b] = t2 + 0.5 * (p - T1[i,b] * T1[j,a])
        end
        
        # Form Intermediates F and W
        # F_ae
        F_ae = zeros(Float64, Nv, Nv)
        for e in 1:Nv, a in 1:Nv
            val = (a==e ? f[No+a] : 0.0)
            # -0.5 sum_m f_me t_m^a ... f_me is 0 for canonical HF
            # sum_mf t_m^f <ma||fe>
            sum_mf = 0.0
            for m in 1:No, f in 1:Nv
                sum_mf += T1[m,f] * get_v(m, No+a, No+f, No+e)
            end
            val += sum_mf
            # -0.5 sum_mef t_ma^ef <ma||ef>
            sum_mef = 0.0
            for m in 1:No, n in 1:No, f in 1:Nv
                 sum_mef += Tau_tilde[m,n,a,f] * get_v(m, n, No+e, No+f)
            end
            val -= 0.5 * sum_mef
            F_ae[a,e] = val
        end
        
        # F_mi
        F_mi = zeros(Float64, No, No)
        for i in 1:No, m in 1:No
            val = (m==i ? f[m] : 0.0)
            # + sum_ne t_n^e <mn||ie>
            sum_ne = 0.0
            for n in 1:No, e in 1:Nv
                sum_ne += T1[n,e] * get_v(m, n, i, No+e)
            end
            val += sum_ne
            # + 0.5 sum_nef t_in^ef <mn||ef>
            sum_nef = 0.0
            for n in 1:No, e in 1:Nv, f in 1:Nv
                sum_nef += Tau_tilde[i,n,e,f] * get_v(m, n, No+e, No+f)
            end
            val += 0.5 * sum_nef
            F_mi[m,i] = val
        end
        
        # F_me
        F_me = zeros(Float64, No, Nv)
        for e in 1:Nv, m in 1:No
            val = 0.0 # f_me is 0
            # sum_nf t_n^f <mn||ef>
            for n in 1:No, f in 1:Nv
                val += T1[n,f] * get_v(m, n, No+e, No+f)
            end
            F_me[m,e] = val
        end
        
        # W_mnij
        W_mnij = zeros(Float64, No, No, No, No)
        for j in 1:No, i in 1:No, n in 1:No, m in 1:No
            val = get_v(m, n, i, j)
            # P(ij) sum_e t_j^e <mn||ie>
            val += sum(T1[j,e] * get_v(m, n, i, No+e) for e in 1:Nv)
            val -= sum(T1[i,e] * get_v(m, n, j, No+e) for e in 1:Nv)
            # + 0.25 sum_ef tau_ij^ef <mn||ef>
            val += 0.25 * sum(Tau[i,j,e,f] * get_v(m, n, No+e, No+f) for e in 1:Nv, f in 1:Nv)
            W_mnij[m,n,i,j] = val
        end
        
        # W_abef
        W_abef = zeros(Float64, Nv, Nv, Nv, Nv)
        for f in 1:Nv, e in 1:Nv, b in 1:Nv, a in 1:Nv
            val = get_v(No+a, No+b, No+e, No+f)
            # - P(ab) sum_m t_m^b <ma||ef>
            val -= sum(T1[m,b] * get_v(m, No+a, No+e, No+f) for m in 1:No)
            val += sum(T1[m,a] * get_v(m, No+b, No+e, No+f) for m in 1:No)
            # + 0.25 sum_mn tau_mn^ab <mn||ef>
            val += 0.25 * sum(Tau[m,n,a,b] * get_v(m, n, No+e, No+f) for m in 1:No, n in 1:No)
            W_abef[a,b,e,f] = val
        end
        
        # W_mbej
        W_mbej = zeros(Float64, No, Nv, Nv, No)
        for j in 1:No, e in 1:Nv, b in 1:Nv, m in 1:No
            val = get_v(m, No+b, No+e, j)
            # + sum_f t_j^f <mb||ef>
            val += sum(T1[j,f] * get_v(m, No+b, No+e, No+f) for f in 1:Nv)
            # - sum_n t_n^b <mn||ej>
            val -= sum(T1[n,b] * get_v(m, n, No+e, j) for n in 1:No)
            # - 0.5 sum_nf (t_jn^fb + 2 t_j^f t_n^b) <mn||ef>
            # Standard: -0.5 sum_nf T_jn^fb <mn||ef> - t_j^f t_n^b <mn||ef> ?
            # Using Tau: -0.5 sum_nf (t_jn^fb + 2 t_j^f t_n^b)
            # Actually easier: - sum_nf (0.5 * T2[j,n,f,b] + T1[j,f]*T1[n,b]) * <mn||ef>
            val -= sum( (0.5 * T2[j,n,f,b] + T1[j,f]*T1[n,b]) * get_v(m, n, No+e, No+f) for n in 1:No, f in 1:Nv)
            W_mbej[m,b,e,j] = val
        end
        
        # Residuals T1
        # r_i^a
        R1 = zeros(Float64, No, Nv)
        for a in 1:Nv, i in 1:No
            # f_ia (0)
            val = 0.0 # f[i, No+a]
            # + sum_e t_i^e F_ae
            val += sum(T1[i,e] * F_ae[a,e] for e in 1:Nv)
            # - sum_m t_m^a F_mi
            val -= sum(T1[m,a] * F_mi[m,i] for m in 1:No)
            # + sum_me t_im^ae F_me
            val += sum(T2[i,m,a,e] * F_me[m,e] for m in 1:No, e in 1:Nv)
            # - sum_nf t_n^f <na||if>
            val -= sum(T1[n,f] * get_v(n, No+a, i, No+f) for n in 1:No, f in 1:Nv)
            # - 0.5 sum_mef t_im^ef <ma||ef>
            val -= 0.5 * sum(T2[i,m,e,f] * get_v(m, No+a, No+e, No+f) for m in 1:No, e in 1:Nv, f in 1:Nv)
            # - 0.5 sum_mne t_mn^ae <nm||ie>
            val -= 0.5 * sum(T2[m,n,a,e] * get_v(n, m, i, No+e) for m in 1:No, n in 1:No, e in 1:Nv)
            R1[i,a] = val
        end
        
        # Residuals T2
        # r_ij^ab
        R2 = zeros(Float64, No, No, Nv, Nv)
        for b in 1:Nv, a in 1:Nv, j in 1:No, i in 1:No
            # <ij||ab>
            val = get_v(i, j, No+a, No+b)
            # P(ab) sum_e t_ij^ae (F_be - 0.5 sum_m t_m^b F_me)
            term1 = 0.0
            for e in 1:Nv
                term1 += T2[i,j,a,e] * (F_ae[b,e] - 0.5 * sum(T1[m,b] * F_me[m,e] for m in 1:No))
            end
            # Permute a,b
            term1_perm = 0.0
            for e in 1:Nv
                term1_perm += T2[i,j,b,e] * (F_ae[a,e] - 0.5 * sum(T1[m,a] * F_me[m,e] for m in 1:No))
            end
            val += term1 - term1_perm
            
            # - P(ij) sum_m t_im^ab (F_mj + 0.5 sum_e t_j^e F_me)
            term2 = 0.0
            for m in 1:No
                term2 += T2[i,m,a,b] * (F_mi[m,j] + 0.5 * sum(T1[j,e] * F_me[m,e] for e in 1:Nv))
            end
            term2_perm = 0.0
            for m in 1:No
                term2_perm += T2[j,m,a,b] * (F_mi[m,i] + 0.5 * sum(T1[i,e] * F_me[m,e] for e in 1:Nv))
            end
            val -= (term2 - term2_perm)
            
            # + 0.5 sum_mn tau_mn^ab W_mnij
            val += 0.5 * sum(Tau[m,n,a,b] * W_mnij[m,n,i,j] for m in 1:No, n in 1:No)
            
            # + 0.5 sum_ef tau_ij^ef W_abef
            val += 0.5 * sum(Tau[i,j,e,f] * W_abef[a,b,e,f] for e in 1:Nv, f in 1:Nv)
            
            # + P(ij) P(ab) sum_me (t_im^ae W_mbej - t_i^e t_m^a <mb||ej>)  <-- Correction needed for W definition
            # Let's expand P(ij)P(ab) = (1 - P_ij)(1 - P_ab) = 1 - P_ab - P_ij + P_ij P_ab
            # Loop for 4 permutations
            
            # Base term(i,j,a,b)
            t_base = sum( (T2[i,m,a,e] * W_mbej[m,b,e,j] - T1[i,e]*T1[m,a]*get_v(m, No+b, No+e, j)) for m in 1:No, e in 1:Nv)
            # Perm a,b
            t_p_ab = sum( (T2[i,m,b,e] * W_mbej[m,a,e,j] - T1[i,e]*T1[m,b]*get_v(m, No+a, No+e, j)) for m in 1:No, e in 1:Nv)
            # Perm i,j
            t_p_ij = sum( (T2[j,m,a,e] * W_mbej[m,b,e,i] - T1[j,e]*T1[m,a]*get_v(m, No+b, No+e, i)) for m in 1:No, e in 1:Nv)
            # Perm both
            t_p_all = sum( (T2[j,m,b,e] * W_mbej[m,a,e,i] - T1[j,e]*T1[m,b]*get_v(m, No+a, No+e, i)) for m in 1:No, e in 1:Nv)
            
            val += (t_base - t_p_ab - t_p_ij + t_p_all)
            
            # + P(ij) sum_e t_i^e <ab||ej>
            val += sum(T1[i,e] * get_v(No+a, No+b, No+e, j) for e in 1:Nv)
            val -= sum(T1[j,e] * get_v(No+a, No+b, No+e, i) for e in 1:Nv)
            
            # - P(ab) sum_m t_m^a <mb||ij>
            val -= sum(T1[m,a] * get_v(m, No+b, i, j) for m in 1:No)
            val += sum(T1[m,b] * get_v(m, No+a, i, j) for m in 1:No)
            
            R2[i,j,a,b] = val
        end
        
        # Calculate Energy
        # E_corr = sum_ia f_ia t_i^a + 0.25 sum_ijab <ij||ab> t_ij^ab + 0.5 sum_ijab <ij||ab> t_i^a t_j^b
        E_corr = 0.0
        # f_ia is 0 in canonical HF
        for i in 1:No, j in 1:No, a in 1:Nv, b in 1:Nv
            v_int = get_v(i, j, No+a, No+b)
            E_corr += 0.25 * v_int * T2[i,j,a,b] + 0.5 * v_int * T1[i,a] * T1[j,b]
        end
        
        delta_E = E_corr - E_corr_old
        E_corr_old = E_corr
        
        # Update Amplitudes
        T1_new = T1 .+ R1 ./ D1
        T2_new = T2 .+ R2 ./ D2
        
        # DIIS
        # DIISManager expects vectors.
        # But we have T1 and T2. We can pass a tuple.
        # Errors are R1 and R2.
        
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

function RunRCCSD(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int)
    println("--- Starting RCCSD Calculation ---")
    println("Running RHF...")
    RHF_Res = ChouChem.RunRHF(MolInAng, Charge, Multiplicity)
    
    # Convert to UHF structure (Spin Orbitals)
    # This just doubles the basis and data.
    println("Converting RHF to Spin-Orbitals...")
    UHF_Res = ChouChem.RHF2UHF(RHF_Res)
    
    # Run CCSD
    return spin_orbital_ccsd(UHF_Res)
end
