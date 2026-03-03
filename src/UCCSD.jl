export RunUCCSD

function RunUCCSD(Molecule::Vector{Atom}, Charge::Int, Multiplicity::Int; MaxIter = 100, Threshold = 1e-7, r_conv = 1e-7, max_diis = 8)
    uhf = RunUHF(Molecule, Charge, Multiplicity; MaxIter = MaxIter, Threshold = Threshold)
    return RunUCCSD(uhf; maxiter = MaxIter, e_conv = Threshold, r_conv = r_conv, max_diis = max_diis)
end

function RunUCCSD(uhf::UHFResults; maxiter = 100, e_conv = 1e-7, r_conv = 1e-7, max_diis = 8)
    println("--- Starting GSO-UCCSD Calculation ---")

    # 1. Integral Transformation and GSO Construction
    # ==========================================================
    norb = size(uhf.CAlpha, 1)
    nocc_a = uhf.ENumAlpha
    nocc_b = uhf.ENumBeta
    nvir_a = norb - nocc_a
    nvir_b = norb - nocc_b
    
    nocc_tot = nocc_a + nocc_b
    nvir_tot = nvir_a + nvir_b
    ntot = 2 * norb

    println("Transforming integrals to GSO basis...")
    
    # Transform to MO basis (Chemist notation <pq|rs>)
    mo_aaaa = TransERI_blas(uhf.ERI, uhf.CAlpha, uhf.CAlpha, uhf.CAlpha, uhf.CAlpha, norb)
    mo_bbbb = TransERI_blas(uhf.ERI, uhf.CBeta, uhf.CBeta, uhf.CBeta, uhf.CBeta, norb)
    mo_aabb = TransERI_blas(uhf.ERI, uhf.CAlpha, uhf.CAlpha, uhf.CBeta, uhf.CBeta, norb)
    mo_bbaa = permutedims(mo_aabb, (3, 4, 1, 2))

    # Create permutation map to Occ-Vir order
    # Target: [OccA, OccB, VirA, VirB]
    perm_map = zeros(Int, ntot)
    for i in 1:nocc_a; perm_map[i] = i; end
    for i in 1:nocc_b; perm_map[nocc_a + i] = norb + i; end
    for i in 1:nvir_a; perm_map[nocc_tot + i] = nocc_a + i; end
    for i in 1:nvir_b; perm_map[nocc_tot + nvir_a + i] = norb + nocc_b + i; end

    # Build full GSO tensor in raw block order
    ERI_raw = zeros(Float64, ntot, ntot, ntot, ntot)
    ERI_raw[1:norb, 1:norb, 1:norb, 1:norb] .= mo_aaaa
    ERI_raw[norb+1:2norb, norb+1:2norb, norb+1:2norb, norb+1:2norb] .= mo_bbbb
    ERI_raw[1:norb, 1:norb, norb+1:2norb, norb+1:2norb] .= mo_aabb
    ERI_raw[norb+1:2norb, norb+1:2norb, 1:norb, 1:norb] .= mo_bbaa

    # Permute to Occ-Vir order <pq|rs> Chemist
    ERI_chem = ERI_raw[perm_map, perm_map, perm_map, perm_map]
    
    # Convert to Physicist notation <pq||rs> = <pq|rs> - <pq|sr>
    # Phys <pq|rs> = Chem <pr|qs>
    ERI_phys_raw = permutedims(ERI_chem, (1, 3, 2, 4))
    ERI = ERI_phys_raw .- permutedims(ERI_phys_raw, (1, 2, 4, 3))
    
    # Free memory
    mo_aaaa = nothing; mo_bbbb = nothing; mo_aabb = nothing; mo_bbaa = nothing
    ERI_raw = nothing; ERI_chem = nothing; ERI_phys_raw = nothing
    GC.gc()

    # Fock Matrix (Diagonal energies in canonical HF)
    E_raw = zeros(Float64, ntot)
    E_raw[1:norb] .= uhf.EAlpha
    E_raw[norb+1:2norb] .= uhf.EBeta
    E_sorted = E_raw[perm_map]
    
    o = 1:nocc_tot
    v = nocc_tot+1:ntot
    fock_diag = E_sorted
    
    # Precompute Denominators
    eps_o = E_sorted[o]
    eps_v = E_sorted[v]
    D1 = [eps_o[i] - eps_v[a] for i in 1:nocc_tot, a in 1:nvir_tot]
    D2 = [eps_o[i] + eps_o[j] - eps_v[a] - eps_v[b] for i in 1:nocc_tot, j in 1:nocc_tot, a in 1:nvir_tot, b in 1:nvir_tot]

    # 2. Initial Amplitudes (MP2 Guess)
    # ==========================================================
    t1 = zeros(Float64, nocc_tot, nvir_tot)
    t2 = zeros(Float64, nocc_tot, nocc_tot, nvir_tot, nvir_tot)
    
    oovv = ERI[o, o, v, v]
    t2 .= oovv ./ D2
    
    E_mp2 = 0.25 * dot(t2, oovv)
    println(@sprintf("Init T2, MP2 Energy = %.15f", E_mp2))

    Ecc = E_mp2
    
    # Extract blocks for efficiency
    oooo = ERI[o, o, o, o]
    ooov = ERI[o, o, o, v]
    # oovv is already defined
    ovov = ERI[o, v, o, v]
    ovvv = ERI[o, v, v, v]
    vvvv = ERI[v, v, v, v]
    
    diis = DIISManager{Tuple{Matrix{Float64}, Array{Float64, 4}}}(max_diis)
    
    # 3. Iteration Loop
    # ==========================================================
    for iter in 1:maxiter
        Ecc_last = Ecc
        
        t1new, t2new = update_amps(t1, t2, fock_diag, o, v, oooo, ooov, oovv, ovov, ovvv, vvvv, D1, D2)
        
        # Calculate Energy
        # E = 0.25 * sum t2_ijab * V_ijab + 0.5 * sum t1_ia * t1_jb * V_ijab
        E_t2 = 0.25 * dot(t2, oovv)
        
        vec_t1 = reshape(t1, nocc_tot * nvir_tot)
        mat_oovv = reshape(permutedims(oovv, (1,3,2,4)), (nocc_tot*nvir_tot, nocc_tot*nvir_tot))
        E_t1 = 0.5 * dot(vec_t1, mat_oovv * vec_t1)
        
        Ecc = E_t2 + E_t1
        dE = Ecc - Ecc_last
        
        diff_t1 = t1new .- t1
        diff_t2 = t2new .- t2
        rms = sqrt(dot(diff_t1, diff_t1) + dot(diff_t2, diff_t2))
        
        println(@sprintf("CC Iter %3d: CC Ecorr = %.15f  dE = % .5E  rms = % .5E", iter, Ecc, dE, rms))
        
        if abs(dE) < e_conv && abs(rms) < r_conv
            println("\nGSO-UCCSD converged.")
            return CCSDResults((uhf.ENumAlpha + uhf.ENumBeta), nocc_tot, nvir_tot, t1new, t2new, Ecc, uhf.Etot + Ecc)
        end
        
        t1, t2 = diis_update!(diis, (t1new, t2new), (diff_t1, diff_t2))
    end
    
    println("GSO-UCCSD failed to converge.")
    return CCSDResults((uhf.ENumAlpha + uhf.ENumBeta), nocc_tot, nvir_tot, t1, t2, Ecc, uhf.Etot + Ecc)
end

function update_amps(t1, t2, fock_diag, o_idx, v_idx, oooo, ooov, oovv, ovov, ovvv, vvvv, D1, D2)
    nocc, nvir = size(t1)
    ev = fock_diag[v_idx]
    eo = fock_diag[o_idx]
    
    # 1. Precompute Tau
    # tau_ijab = t2_ijab + 0.5 * (t1_ia*t1_jb - t1_ib*t1_ja)
    tau = copy(t2)
    for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc
        tau[i,j,a,b] += 0.5 * (t1[i,a]*t1[j,b] - t1[i,b]*t1[j,a])
    end
    # tau_tilde_ijab = t2_ijab + 0.25 * (t1_ia*t1_jb - t1_ib*t1_ja)
    tau_tilde = copy(t2)
    for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc
        tau_tilde[i,j,a,b] += 0.25 * (t1[i,a]*t1[j,b] - t1[i,b]*t1[j,a])
    end
    
    # 2. Intermediates F
    # Fov_me = sum_nf t1_nf * <mn||ef>
    # oovv is (m,n,e,f). Contract n,f.
    mat_oovv_nf = reshape(permutedims(oovv, (1,3,2,4)), (nocc*nvir, nocc*nvir)) # (me, nf)
    vec_t1 = reshape(t1, nocc*nvir) # (nf)
    Fov = reshape(mat_oovv_nf * vec_t1, (nocc, nvir)) # (me)
    
    # Fvv_ae = f_ae - 0.5 sum_m t1_ma Fov_me + sum_mf t1_mf <am||ef> - 0.5 sum_mnf tau_tilde_mnaf <mn||ef>
    Fvv = diagm(ev)
    Fvv .-= 0.5 .* (t1' * Fov) # sum_m t1_ma * Fov_me -> (a, e)
    # Term 2: sum_mf t1_mf <am||ef> = - sum_mf t1_mf <ma||ef> = - sum_mf t1_mf * ovvv[m,a,e,f]
    # ovvv is (m,a,e,f). Contract m,f.
    mat_ovvv_mf = reshape(permutedims(ovvv, (2,3,1,4)), (nvir*nvir, nocc*nvir)) # (ae, mf)
    Fvv .-= reshape(mat_ovvv_mf * vec_t1, (nvir, nvir))
    # Term 3: -0.5 sum_mnf tau_tilde_mnaf * oovv_mnef
    # tau_tilde is (m,n,a,f). oovv is (m,n,e,f). Contract m,n,f.
    mat_tau_mnf = reshape(permutedims(tau_tilde, (3, 1, 2, 4)), (nvir, nocc*nocc*nvir)) # (a, mnf)
    mat_oovv_mnf = reshape(permutedims(oovv, (3, 1, 2, 4)), (nvir, nocc*nocc*nvir)) # (e, mnf)
    Fvv .-= 0.5 .* (mat_tau_mnf * mat_oovv_mnf')
    
    # Foo_mi = f_mi + 0.5 sum_e t1_ie Fov_me + sum_ne t1_ne <mn||ie> + 0.5 sum_nef tau_tilde_inef <mn||ef>
    Foo = diagm(eo)
    Foo .+= 0.5 .* (Fov * t1') # sum_e Fov_me * t1_ie -> (m, i)
    # Term 2: sum_ne t1_ne * ooov_mnie
    # ooov is (m,n,i,e). Contract n,e.
    mat_ooov_ne = reshape(permutedims(ooov, (1,3,2,4)), (nocc*nocc, nocc*nvir)) # (mi, ne)
    Foo .+= reshape(mat_ooov_ne * vec_t1, (nocc, nocc))
    # Term 3: 0.5 sum_nef tau_tilde_inef * oovv_mnef
    # tau_tilde is (i,n,e,f). oovv is (m,n,e,f). Contract n,e,f.
    mat_tau_nef = reshape(tau_tilde, (nocc, nocc*nvir*nvir)) # (i, nef)
    mat_oovv_nef = reshape(oovv, (nocc, nocc*nvir*nvir)) # (m, nef)
    Foo .+= 0.5 .* (mat_oovv_nef * mat_tau_nef') # (m, i)
    
    # 3. Intermediates W
    # Woooo_mnij = <mn||ij> + P(ij) sum_e t1_je <mn||ie> + 0.25 sum_ef tau_ijef <mn||ef>
    Woooo = copy(oooo)
    # Term 2: sum_e (t1_je * ooov_mnie - t1_ie * ooov_mnje)
    mat_ooov_e = reshape(ooov, (nocc*nocc*nocc, nvir)) # (mni, e)
    term2_oooo = reshape(mat_ooov_e * t1', (nocc, nocc, nocc, nocc)) # (mni, j)
    Woooo .+= term2_oooo .- permutedims(term2_oooo, (1, 2, 4, 3))
    # Term 3: 0.25 sum_ef oovv_mnef * tau_ijef
    mat_tau_ef = reshape(tau, (nocc*nocc, nvir*nvir)) # (ij, ef)
    mat_oovv_ef = reshape(oovv, (nocc*nocc, nvir*nvir)) # (mn, ef)
    Woooo .+= 0.25 .* reshape(mat_oovv_ef * mat_tau_ef', (nocc, nocc, nocc, nocc))
    
    # Wvvvv_abef = <ab||ef> - P(ab) sum_m t1_mb <am||ef> + 0.25 sum_mn tau_mnab <mn||ef>
    Wvvvv = copy(vvvv)
    # Term 2: sum_m (t1_mb * ovvv_maef - t1_ma * ovvv_mbef)
    # Note: <am||ef> = -ovvv[m,a,e,f]. So term is + sum_m t1_mb * ovvv_maef
    mat_ovvv_m = reshape(ovvv, (nocc, nvir*nvir*nvir)) # (m, aef)
    term2_vvvv = reshape(t1' * mat_ovvv_m, (nvir, nvir, nvir, nvir)) # (b, a, e, f)
    Wvvvv .+= permutedims(term2_vvvv, (2, 1, 3, 4)) .- permutedims(term2_vvvv, (1, 2, 3, 4))
    # Term 3: 0.25 sum_mn tau_mnab * oovv_mnef
    Wvvvv .+= 0.25 .* reshape(mat_tau_ef' * mat_oovv_ef, (nvir, nvir, nvir, nvir))
    
    # Wovvo_mbej = <mb||ej> + sum_f t1_jf <mb||ef> - sum_n t1_nb <mn||ej> - (0.5 sum_nf t2_jnfb + sum_f t1_jf sum_n t1_nb) <mn||ef>
    # <mb||ej> = -<mb||je> = -ovov[m,b,j,e]. Permute to (m,b,e,j).
    Wovvo = -1.0 .* permutedims(ovov, (1, 2, 4, 3))
    # Term 2: sum_f t1_jf * ovvv_mbef
    mat_ovvv_f = reshape(ovvv, (nocc*nvir*nvir, nvir)) # (mbe, f)
    Wovvo .+= reshape(mat_ovvv_f * t1', (nocc, nvir, nvir, nocc)) # (mbe, j)
    # Term 3: - sum_n t1_nb * <mn||ej> = sum_n t1_nb * <mn||je> = sum_n t1_nb * ooov_mnje
    mat_ooov_n = reshape(permutedims(ooov, (1, 3, 4, 2)), (nocc*nocc*nvir, nocc)) # (mje, n)
    term3_Wovvo = reshape(mat_ooov_n * t1, (nocc, nocc, nvir, nvir)) # (mje, b) -> (m, j, e, b)
    Wovvo .+= permutedims(term3_Wovvo, (1, 4, 3, 2)) # (m, b, e, j)
    # Term 4: -0.5 sum_nf t2_jnfb * oovv_mnef
    # Contract n,f. t2(j,n,f,b) -> (j,b,n,f). oovv(m,n,e,f) -> (m,e,n,f).
    mat_t2_nf = reshape(permutedims(t2, (1, 4, 2, 3)), (nocc*nvir, nocc*nvir)) # (jb, nf)
    Wovvo .-= 0.5 .* permutedims(reshape(mat_t2_nf * mat_oovv_nf', (nocc, nvir, nocc, nvir)), (3, 2, 4, 1))
    # Term 5: - sum_f t1_jf * (sum_n t1_nb * oovv_mnef)
    mat_oovv_n = reshape(permutedims(oovv, (1, 3, 4, 2)), (nocc*nvir*nvir, nocc)) # (mef, n)
    inner_Wovvo = reshape(mat_oovv_n * t1, (nocc, nvir, nvir, nvir)) # (mef, b) -> (m, e, f, b)
    mat_inner_W = reshape(permutedims(inner_Wovvo, (1, 2, 4, 3)), (nocc*nvir*nvir, nvir)) # (meb, f)
    res_Wovvo = reshape(mat_inner_W * t1', (nocc, nvir, nvir, nocc)) # (meb, j) -> (m, e, b, j)
    Wovvo .-= permutedims(res_Wovvo, (1, 3, 2, 4)) # (m, b, e, j)
    
    # 4. Update Residuals
    # Singles (T1) Residual
    # r_ia = f_ia + sum_e t1_ie Fvv_ae - sum_m t1_ma Foo_mi + sum_me t2_imae Fov_me - sum_nf t1_nf <na||if> - 0.5 sum_mef t2_imef <ma||ef> - 0.5 sum_mne t2_mnae <mn||ie>
    t1new = t1 * Fvv' .- Foo' * t1
    # Term 3: sum_me t2_imae Fov_me
    mat_t2_me = reshape(permutedims(t2, (1, 3, 2, 4)), (nocc*nvir, nocc*nvir)) # (ia, me)
    t1new .+= reshape(mat_t2_me * reshape(Fov, nocc*nvir), (nocc, nvir))
    # Term 4: - sum_nf t1_nf * <na||if> = sum_nf t1_nf * <na||fi> = sum_nf t1_nf * ovov_nafi
    mat_ovov_nf = reshape(permutedims(ovov, (2, 3, 1, 4)), (nvir*nocc, nocc*nvir)) # (ai, nf)
    t1new .+= reshape(mat_ovov_nf * vec_t1, (nvir, nocc))' # (ai -> ia)
    # Term 5: -0.5 sum_mef t2_imef * <ma||ef> = -0.5 sum_mef t2_imef * ovvv_maef
    mat_t2_mef = reshape(t2, (nocc, nocc*nvir*nvir)) # (i, mef)
    mat_ovvv_mef = reshape(permutedims(ovvv, (2, 1, 3, 4)), (nvir, nocc*nvir*nvir)) # (a, mef)
    t1new .-= 0.5 .* (mat_t2_mef * mat_ovvv_mef')
    # Term 6: -0.5 sum_mne t2_mnae * <mn||ie> = -0.5 sum_mne t2_mnae * ooov_mnie
    mat_t2_mne = reshape(permutedims(t2, (3, 1, 2, 4)), (nvir, nocc*nocc*nvir)) # (a, mne)
    mat_ooov_mne = reshape(permutedims(ooov, (3, 1, 2, 4)), (nocc, nocc*nocc*nvir)) # (i, mne)
    t1new .-= 0.5 .* (mat_ooov_mne * mat_t2_mne')
    
    # Doubles (T2) Residual
    # r_ijab = <ij||ab> + P(ab) sum_e t2_ijae Fvv_tmp_be - P(ij) sum_m t2_imab Foo_tmp_mj + 0.5 sum_mn tau_mnab Woooo_mnij + 0.5 sum_ef tau_ijef Wvvvv_abef + P(ij)P(ab) (sum_me t2_imae Wovvo_mbej - sum_me t1_ie t1_ma ovov_mbje) + P(ab) sum_e t1_ie <ab||ej> - P(ij) sum_m t1_ma <mb||ij>
    t2new = copy(oovv)
    Fvv_tmp = Fvv .- 0.5 .* (t1' * Fov)
    Foo_tmp = Foo .+ 0.5 .* (Fov * t1')
    # Term 2: P(ab) sum_e t2_ijae * Fvv_tmp_be
    mat_t2_e = reshape(t2, (nocc*nocc*nvir, nvir)) # (ija, e)
    term_vv = reshape(mat_t2_e * Fvv_tmp', (nocc, nocc, nvir, nvir)) # (ija, b)
    t2new .+= term_vv .- permutedims(term_vv, (1, 2, 4, 3))
    # Term 3: - P(ij) sum_m t2_imab * Foo_tmp_mj
    mat_t2_m = reshape(permutedims(t2, (1, 3, 4, 2)), (nocc*nvir*nvir, nocc)) # (iab, m)
    term_oo = reshape(mat_t2_m * Foo_tmp, (nocc, nvir, nvir, nocc)) # (iab, j)
    term_oo_perm = permutedims(term_oo, (1, 4, 2, 3)) # (i, j, a, b)
    t2new .-= term_oo_perm .- permutedims(term_oo_perm, (2, 1, 3, 4))
    # Term 4: 0.5 sum_mn tau_mnab * Woooo_mnij
    t2new .+= 0.5 .* permutedims(reshape(mat_tau_ef' * reshape(Woooo, (nocc*nocc, nocc*nocc)), (nvir, nvir, nocc, nocc)), (3, 4, 1, 2))
    # Term 5: 0.5 sum_ef tau_ijef * Wvvvv_abef
    t2new .+= 0.5 .* reshape(mat_tau_ef * reshape(Wvvvv, (nvir*nvir, nvir*nvir))', (nocc, nocc, nvir, nvir))
    # Term 6: P(ij)P(ab) (sum_me t2_imae Wovvo_mbej - sum_me t1_ie t1_ma ovov_mbje)
    # Note: <mb||je> = ovov[m,b,j,e].
    # Part A: sum_me t2_imae * Wovvo_mbej
    mat_t2_me_A = reshape(permutedims(t2, (1, 3, 2, 4)), (nocc*nvir, nocc*nvir)) # (ia, me)
    mat_Wovvo_me = reshape(permutedims(Wovvo, (2, 4, 1, 3)), (nvir*nocc, nocc*nvir)) # (bj, me)
    term_A = reshape(mat_t2_me_A * mat_Wovvo_me', (nocc, nvir, nvir, nocc)) # (ia, bj) -> (i, a, b, j)
    # Part B: sum_me t1_ie * t1_ma * ovov_mbje
    mat_ovov_me = reshape(ovov, (nocc, nvir*nocc*nvir)) # (m, bje)
    inner_B = reshape(t1' * mat_ovov_me, (nvir, nvir, nocc, nvir)) # (a, b, j, e)
    mat_inner_B = reshape(permutedims(inner_B, (1, 2, 3, 4)), (nvir*nvir*nocc, nvir)) # (abj, e)
    res_B = reshape(mat_inner_B * t1', (nvir, nvir, nocc, nocc)) # (abj, i) -> (a, b, j, i)
    # Sign: PySCF is + Part B. (Wait, check gccsd.py again: tmp -= -einsum... yes, +)
    term_mixed_total = permutedims(term_A, (1, 4, 2, 3)) .+ permutedims(res_B, (4, 3, 1, 2))
    # Permute P(ij)P(ab)
    term_mixed_total .-= permutedims(term_mixed_total, (2, 1, 3, 4))
    term_mixed_total .-= permutedims(term_mixed_total, (1, 2, 4, 3))
    t2new .+= term_mixed_total
    # Term 7: P(ab) sum_e t1_ie * <ab||ej> = P(ab) sum_e t1_ie * ovvv_jeba
    mat_ovvv_e_T7 = reshape(permutedims(ovvv, (2, 1, 4, 3)), (nvir, nocc*nvir*nvir)) # (e, jba)
    res_T7 = reshape(t1 * mat_ovvv_e_T7, (nocc, nocc, nvir, nvir)) # (i, jba) -> (i, j, b, a)
    t2new .+= permutedims(res_T7, (1, 2, 4, 3)) .- permutedims(res_T7, (1, 2, 3, 4)) # res_ijab - res_ijba
    # Term 8: - P(ij) sum_m t1_ma * <mb||ij> = - P(ij) sum_m t1_ma * ooov_ijmb
    # Contract m: sum_m ooov_ijmb * t1_ma -> (ijb, a)
    tmp_T8 = permutedims(ooov, (1, 2, 4, 3)) # (i, j, b, m)
    mat_T8 = reshape(tmp_T8, (nocc*nocc*nvir, nocc)) # (ijb, m)
    res_T8 = reshape(mat_T8 * t1, (nocc, nocc, nvir, nvir)) # (ijb, a) -> (i, j, b, a)
    t2new .-= permutedims(res_T8, (1, 2, 4, 3)) .- permutedims(res_T8, (2, 1, 4, 3)) # res_ijab - res_jiab
    
    t1new ./= D1
    t2new ./= D2
    return t1new, t2new
end