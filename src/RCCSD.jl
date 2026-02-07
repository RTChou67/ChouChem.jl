export RunRCCSD


function RunRCCSD(Molecule::Vector{Atom}, Charge::Int, Multiplicity::Int; MaxIter = 100, Threshold = 1e-7, r_conv = 1e-7, max_diis = 8)
	rhf = RunRHF(Molecule, Charge, Multiplicity; MaxIter = MaxIter, Threshold = Threshold)
	return RunRCCSD(rhf; maxiter = MaxIter, e_conv = Threshold, r_conv = r_conv, max_diis = max_diis)
end

function RunRCCSD(rhf::RHFResults; maxiter = 100, e_conv = 1e-7, r_conv = 1e-7, max_diis = 8)
	println("--- Starting RCCSD Calculation ---")

	norb = size(rhf.C, 1)
	nocc = rhf.ENum ÷ 2
	nvir = norb - nocc

	println("MO Transformation...")
	C = rhf.C
	ERI_MO = TransERI_blas(rhf.ERI, C, C, C, C, norb)

	ERI = permutedims(ERI_MO, (1, 3, 2, 4))
	L = 2.0 .* ERI .- permutedims(ERI, (1, 2, 4, 3))

	o = 1:nocc
	v = (nocc+1):norb

	eps = rhf.E
	F = diagm(eps)

	Dia = [eps[i] - eps[a] for i in o, a in v]
	Dijab = [eps[i] + eps[j] - eps[a] - eps[b] for i in o, j in o, a in v, b in v]

	t1 = zeros(Float64, nocc, nvir)
	t2 = zeros(Float64, nocc, nocc, nvir, nvir)

	ERI_oovv = ERI[o, o, v, v]
	t2 .= ERI_oovv ./ Dijab

	Ecc = cc_energy(o, v, F, L, t1, t2)
	println(@sprintf("CC Iter %3d: CC Ecorr = %.15f  dE = % .5E", 0, Ecc, -Ecc))

	diis = DIISManager{Tuple{Matrix{Float64}, Array{Float64, 4}}}(max_diis)

	for iter in 1:maxiter
		Ecc_last = Ecc
		r1, r2 = residuals(o, v, F, ERI, L, t1, t2)

		t1_new = t1 .+ r1 ./ Dia
		t2_new = t2 .+ r2 ./ Dijab

		diff_t1 = (t1_new .- t1)
		diff_t2 = (t2_new .- t2)
		rms = sqrt(dot(diff_t1, diff_t1) + dot(diff_t2, diff_t2))

		t1, t2 = t1_new, t2_new
		Ecc = cc_energy(o, v, F, L, t1, t2)
		dE = Ecc - Ecc_last

		println(@sprintf("CC Iter %3d: CC Ecorr = %.15f  dE = % .5E  rms = % .5E", iter, Ecc, dE, rms))

		if abs(dE) < e_conv && abs(rms) < r_conv
			println("\nCCSD converged.")
			return CCSDResults(rhf.ENum, nocc, nvir, t1, t2, Ecc, rhf.Etot + Ecc)
		end
		t1, t2 = diis_update!(diis, (t1, t2), (diff_t1, diff_t2))
	end
	return CCSDResults(rhf.ENum, nocc, nvir, t1, t2, Ecc, rhf.Etot + Ecc)
end

function cc_energy(o, v, F, L, t1, t2)
	E1 = 2.0 * dot(F[o, v], t1)
	L_oovv = L[o, o, v, v]
	E2_a = dot(t2, L_oovv)
	E2_b = 0.0
	nocc, nvir = size(t1)
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc
		E2_b += t1[i, a] * t1[j, b] * L_oovv[i, j, a, b]
	end
	return E1 + E2_a + E2_b
end

function residuals(o, v, F, ERI, L, t1, t2)
	Fae = build_Fae(o, v, F, L, t1, t2)
	Fmi = build_Fmi(o, v, F, L, t1, t2)
	Fme = build_Fme(o, v, F, L, t1)
	Wmnij = build_Wmnij(o, v, ERI, t1, t2)
	Wmbej = build_Wmbej(o, v, ERI, L, t1, t2)
	Wmbje = build_Wmbje(o, v, ERI, t1, t2)
	Zmbij = build_Zmbij(o, v, ERI, t1, t2)
	r1 = r_T1(o, v, F, ERI, L, t1, t2, Fae, Fme, Fmi)
	r2 = r_T2(o, v, F, ERI, L, t1, t2, Fae, Fme, Fmi, Wmnij, Wmbej, Wmbje, Zmbij)
	return r1, r2
end

function build_Fae(o, v, F, L, t1, t2)
	Fae = F[v, v] .- 0.5 .* permutedims((F[o, v]' * t1), (2, 1))
	L_ovvv = L[o, v, v, v]
	nocc, nvir = size(t1)
	temp = zeros(nvir, nvir)
	@inbounds for e in 1:nvir, a in 1:nvir, f in 1:nvir, m in 1:nocc
		temp[a, e] += t1[m, f] * L_ovvv[m, a, f, e]
	end
	Fae .+= temp
	L_oovv = L[o, o, v, v]
	@inbounds for e in 1:nvir, a in 1:nvir, f in 1:nvir, n in 1:nocc, m in 1:nocc
		Fae[a, e] -= (t2[m, n, a, f] + 0.5 * t1[m, a] * t1[n, f]) * L_oovv[m, n, e, f]
	end
	return Fae
end

function build_Fmi(o, v, F, L, t1, t2)
	Fmi = copy(F[o, o])
	Fmi .+= 0.5 .* permutedims((t1 * F[o, v]'), (2, 1))
	L_ooov = L[o, o, o, v]
	nocc, nvir = size(t1)
	temp = zeros(nocc, nocc)
	@inbounds for i in 1:nocc, m in 1:nocc, e in 1:nvir, n in 1:nocc
		temp[m, i] += t1[n, e] * L_ooov[m, n, i, e]
	end
	Fmi .+= temp
	L_oovv = L[o, o, v, v]
	@inbounds for i in 1:nocc, m in 1:nocc, f in 1:nvir, e in 1:nvir, n in 1:nocc
		Fmi[m, i] += (t2[i, n, e, f] + 0.5 * t1[i, e] * t1[n, f]) * L_oovv[m, n, e, f]
	end
	return Fmi
end

function build_Fme(o, v, F, L, t1)
	Fme = copy(F[o, v])
	L_oovv = L[o, o, v, v]
	nocc, nvir = size(t1)
	@inbounds for e in 1:nvir, m in 1:nocc, f in 1:nvir, n in 1:nocc
		Fme[m, e] += t1[n, f] * L_oovv[m, n, e, f]
	end
	return Fme
end

function build_Wmnij(o, v, ERI, t1, t2)
	Wmnij = copy(ERI[o, o, o, o])
	nocc, nvir = size(t1)
	ERI_ooov = ERI[o, o, o, v]
	@inbounds for j in 1:nocc, i in 1:nocc, n in 1:nocc, m in 1:nocc, e in 1:nvir
		Wmnij[m, n, i, j] += t1[j, e] * ERI_ooov[m, n, i, e]
	end
	ERI_oovo = ERI[o, o, v, o]
	@inbounds for j in 1:nocc, i in 1:nocc, n in 1:nocc, m in 1:nocc, e in 1:nvir
		Wmnij[m, n, i, j] += t1[i, e] * ERI_oovo[m, n, e, j]
	end
	ERI_oovv = ERI[o, o, v, v]
	@inbounds for j in 1:nocc, i in 1:nocc, n in 1:nocc, m in 1:nocc, f in 1:nvir, e in 1:nvir
		Wmnij[m, n, i, j] += (t2[i, j, e, f] + t1[i, e] * t1[j, f]) * ERI_oovv[m, n, e, f]
	end
	return Wmnij
end

function build_Wmbej(o, v, ERI, L, t1, t2)
	Wmbej = copy(ERI[o, v, v, o])
	nocc, nvir = size(t1)
	ERI_ovvv = ERI[o, v, v, v]
	@inbounds for j in 1:nocc, e in 1:nvir, b in 1:nvir, m in 1:nocc, f in 1:nvir
		Wmbej[m, b, e, j] += t1[j, f] * ERI_ovvv[m, b, e, f]
	end
	ERI_oovo = ERI[o, o, v, o]
	@inbounds for j in 1:nocc, e in 1:nvir, b in 1:nvir, m in 1:nocc, n in 1:nocc
		Wmbej[m, b, e, j] -= t1[n, b] * ERI_oovo[m, n, e, j]
	end
	ERI_oovv = ERI[o, o, v, v]
	L_oovv = L[o, o, v, v]
	@inbounds for j in 1:nocc, e in 1:nvir, b in 1:nvir, m in 1:nocc, f in 1:nvir, n in 1:nocc
		Wmbej[m, b, e, j] -= (0.5 * t2[j, n, f, b] + t1[j, f] * t1[n, b]) * ERI_oovv[m, n, e, f]
		Wmbej[m, b, e, j] += 0.5 * t2[n, j, f, b] * L_oovv[m, n, e, f]
	end
	return Wmbej
end

function build_Wmbje(o, v, ERI, t1, t2)
	Wmbje = -1.0 .* ERI[o, v, o, v]
	nocc, nvir = size(t1)
	ERI_ovvv = ERI[o, v, v, v]
	@inbounds for e in 1:nvir, j in 1:nocc, b in 1:nvir, m in 1:nocc, f in 1:nvir
		Wmbje[m, b, j, e] -= t1[j, f] * ERI_ovvv[m, b, f, e]
	end
	ERI_ooov = ERI[o, o, o, v]
	@inbounds for e in 1:nvir, j in 1:nocc, b in 1:nvir, m in 1:nocc, n in 1:nocc
		Wmbje[m, b, j, e] += t1[n, b] * ERI_ooov[m, n, j, e]
	end
	ERI_oovv = ERI[o, o, v, v]
	@inbounds for e in 1:nvir, j in 1:nocc, b in 1:nvir, m in 1:nocc, f in 1:nvir, n in 1:nocc
		Wmbje[m, b, j, e] += (0.5 * t2[j, n, f, b] + t1[j, f] * t1[n, b]) * ERI_oovv[m, n, f, e]
	end
	return Wmbje
end

function build_Zmbij(o, v, ERI, t1, t2)
	nocc, nvir = size(t1)
	Zmbij = zeros(nocc, nvir, nocc, nocc)
	ERI_ovvv = ERI[o, v, v, v]
	@inbounds for j in 1:nocc, i in 1:nocc, b in 1:nvir, m in 1:nocc, f in 1:nvir, e in 1:nvir
		Zmbij[m, b, i, j] += ERI_ovvv[m, b, e, f] * (t2[i, j, e, f] + t1[i, e] * t1[j, f])
	end
	return Zmbij
end

function r_T1(o, v, F, ERI, L, t1, t2, Fae, Fme, Fmi)
	r1 = copy(F[o, v])
	nocc, nvir = size(t1)
	r1 .+= t1 * Fae'
	r1 .-= (Fmi' * t1)
	@inbounds for a in 1:nvir, i in 1:nocc, e in 1:nvir, m in 1:nocc
		r1[i, a] += (2.0 * t2[i, m, a, e] - t2[i, m, e, a]) * Fme[m, e]
	end
	L_ovvo = L[o, v, v, o]
	@inbounds for a in 1:nvir, i in 1:nocc, f in 1:nvir, n in 1:nocc
		r1[i, a] += t1[n, f] * L_ovvo[n, a, f, i]
	end
	ERI_ovvv = ERI[o, v, v, v]
	@inbounds for a in 1:nvir, i in 1:nocc, f in 1:nvir, e in 1:nvir, m in 1:nocc
		r1[i, a] += (2.0 * t2[m, i, e, f] - t2[m, i, f, e]) * ERI_ovvv[m, a, e, f]
	end
	L_oovo = L[o, o, v, o]
	@inbounds for a in 1:nvir, i in 1:nocc, e in 1:nvir, n in 1:nocc, m in 1:nocc
		r1[i, a] -= t2[m, n, a, e] * L_oovo[n, m, e, i]
	end
	return r1
end

function r_T2(o, v, F, ERI, L, t1, t2, Fae, Fme, Fmi, Wmnij, Wmbej, Wmbje, Zmbij)
	r2 = 0.5 .* ERI[o, o, v, v]
	nocc, nvir = size(t1)
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, e in 1:nvir
		r2[i, j, a, b] += t2[i, j, a, e] * Fae[b, e]
	end
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, m in 1:nocc
		r2[i, j, a, b] -= t2[i, m, a, b] * Fmi[m, j]
	end
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, n in 1:nocc, m in 1:nocc
		r2[i, j, a, b] += 0.5 * (t2[m, n, a, b] + t1[m, a] * t1[n, b]) * Wmnij[m, n, i, j]
	end
	ERI_vvvv = ERI[v, v, v, v]
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, f in 1:nvir, e in 1:nvir
		r2[i, j, a, b] += 0.5 * (t2[i, j, e, f] + t1[i, e] * t1[j, f]) * ERI_vvvv[a, b, e, f]
	end
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, e in 1:nvir, m in 1:nocc
		r2[i, j, a, b] += (t2[i, m, a, e] - t2[i, m, e, a]) * Wmbej[m, b, e, j]
		r2[i, j, a, b] += t2[i, m, a, e] * (Wmbej[m, b, e, j] + Wmbje[m, b, j, e])
		r2[i, j, a, b] += t2[m, j, a, e] * Wmbje[m, b, i, e]
	end
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, m in 1:nocc
		r2[i, j, a, b] -= t1[m, a] * Zmbij[m, b, i, j]
	end
	ERI_ovvo = ERI[o, v, v, o]
	ERI_ovov = ERI[o, v, o, v]
	ERI_vvvo = ERI[v, v, v, o]
	ERI_ovoo = ERI[o, v, o, o]
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, e in 1:nvir, m in 1:nocc
		r2[i, j, a, b] -= t1[i, e] * t1[m, a] * ERI_ovvo[m, b, e, j]
		r2[i, j, a, b] -= t1[i, e] * t1[m, b] * ERI_ovov[m, a, j, e]
	end
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, e in 1:nvir
		r2[i, j, a, b] += t1[i, e] * ERI_vvvo[a, b, e, j]
	end
	@inbounds for b in 1:nvir, a in 1:nvir, j in 1:nocc, i in 1:nocc, m in 1:nocc
		r2[i, j, a, b] -= t1[m, a] * ERI_ovoo[m, b, i, j]
	end
	r2_final = r2 .+ permutedims(r2, (2, 1, 4, 3))
	return r2_final
end
