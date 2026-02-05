struct CIResults
	EtotHF::Float64
	EtotCI::Float64
	Ecorr::Float64
	StateEnergies::Vector{Float64}
	CIVectors::Matrix{Float64}
	DeterminantSpace::Vector{Vector{Int}}
	MaxExcitation::Int
	NumberOfDeterminants::Int
end



function GenCISpace_Optimized(RefDet::Int, NSpinOrb::Int, MaxExcit::Int)
	println("\n--- Generating CI Space (Optimized) ---")

	Occ = unpack_det(RefDet)
	Vir = setdiff(1:NSpinOrb, Occ)
	nocc = length(Occ)
	nvir = length(Vir)

	Masks = [one(Int) << (i-1) for i in 1:NSpinOrb]

	n_singles = nocc * nvir
	n_doubles = (nocc * (nocc - 1) ÷ 2) * (nvir * (nvir - 1) ÷ 2)

	total_dets = 1
	if MaxExcit >= 1
		;
		total_dets += n_singles;
	end
	if MaxExcit >= 2
		;
		total_dets += n_doubles;
	end

	Dets = Vector{Int}(undef, total_dets)
	Dets[1] = RefDet

	idx = 2

	if MaxExcit >= 1
		@inbounds for i_idx in 1:nocc
			i = Occ[i_idx]
			det_minus_i = RefDet ⊻ Masks[i]
			for a_idx in 1:nvir
				a = Vir[a_idx]
				Dets[idx] = det_minus_i | Masks[a]
				idx += 1
			end
		end
	end

	if MaxExcit >= 2
		@inbounds for i_idx in 1:nocc
			i = Occ[i_idx]
			det_minus_i = RefDet ⊻ Masks[i]

			for j_idx in (i_idx+1):nocc
				j = Occ[j_idx]
				det_minus_ij = det_minus_i ⊻ Masks[j]

				for a_idx in 1:nvir
					a = Vir[a_idx]
					det_plus_a = det_minus_ij | Masks[a]

					for b_idx in (a_idx+1):nvir
						b = Vir[b_idx]
						Dets[idx] = det_plus_a | Masks[b]
						idx += 1
					end
				end
			end
		end
	end

	sort!(Dets)

	@printf("CI Space Size: %d (Ref: 1, S: %d, D: %d)\n", length(Dets), n_singles, n_doubles)
	return Dets
end

function CalcUCI(SCF_Results::UHFResults, MaxExcitation::Int)
	ONum = length(SCF_Results.BasisSet)
	ENum = SCF_Results.ENumAlpha + SCF_Results.ENumBeta
	SONum = 2 * ONum

	C_a, C_b = SCF_Results.CAlpha, SCF_Results.CBeta
	H_a = C_a' * SCF_Results.Hcore * C_a
	H_b = C_b' * SCF_Results.Hcore * C_b

	h_so = zeros(SONum, SONum)
	h_so[1:ONum, 1:ONum] .= H_a
	h_so[(ONum+1):end, (ONum+1):end] .= H_b

	println("\nTransforming ERIs to MO basis...")
	@time begin
		ERI_aa = TransERI_blas(SCF_Results.ERI, C_a, C_a, C_a, C_a, ONum)
		ERI_bb = TransERI_blas(SCF_Results.ERI, C_b, C_b, C_b, C_b, ONum)
		ERI_ab = TransERI_blas(SCF_Results.ERI, C_a, C_a, C_b, C_b, ONum)
		V_so = build_spin_eri(ERI_aa, ERI_bb, ERI_ab, ONum)
	end
	println("ERI transformation complete.")

	RefOcc = vcat(collect(1:SCF_Results.ENumAlpha),
		collect((ONum+1):(ONum+SCF_Results.ENumBeta)))
	RefDet = pack_det(sort(RefOcc))

	if MaxExcitation == -1
		NVir = SONum - ENum
		FCILvl = min(ENum, NVir)
		MaxExcitation = FCILvl
	end

	@printf("Generating CI space for excitation level <= %d...\n", MaxExcitation)
	Dets = GenCISpace_Optimized(RefDet, SONum, MaxExcitation)
	Ndet = length(Dets)

	println("\nBuilding and diagonalizing CI Matrix...")

	I_idx = Int[]
	J_idx = Int[]
	V_val = Float64[]
	sizehint!(I_idx, Ndet * 50)
	sizehint!(J_idx, Ndet * 50)
	sizehint!(V_val, Ndet * 50)

	for i in 1:Ndet
		Di = Dets[i]

		for j in i:Ndet
			Dj = Dets[j]
			diff = Di ⊻ Dj
			n_diff_bits = count_ones(diff)

			if n_diff_bits > 4
				;
				continue;
			end

			val = 0.0

			if n_diff_bits == 0
				occ_i = unpack_det(Di)
				for orb in occ_i
					val += h_so[orb, orb]
				end
				for p_idx in 1:length(occ_i)
					p = occ_i[p_idx]
					for q_idx in (p_idx+1):length(occ_i)
						q = occ_i[q_idx]
						val += V_so[p, p, q, q] - V_so[p, q, q, p]
					end
				end

			elseif n_diff_bits == 2
				p = trailing_zeros(Di & diff) + 1
				q = trailing_zeros(Dj & diff) + 1

				phase = get_phase_single(Di, Dj, p, q)

				val = h_so[p, q]

				common = Di & Dj
				temp_comm = common
				while temp_comm > 0
					k = trailing_zeros(temp_comm) + 1
					val += V_so[p, q, k, k] - V_so[p, k, k, q]
					temp_comm &= ~(one(Int) << (k-1))
				end
				val *= phase

			elseif n_diff_bits == 4
				holes = Di & diff
				p = trailing_zeros(holes) + 1
				q = trailing_zeros(holes & ~(one(Int) << (p-1))) + 1

				parts = Dj & diff
				r = trailing_zeros(parts) + 1
				s = trailing_zeros(parts & ~(one(Int) << (r-1))) + 1

				phase = get_phase_double(Di, Dj, p, q, r, s)

				val = (V_so[p, r, q, s] - V_so[p, s, q, r]) * phase
			end

			if abs(val) > 1e-12
				push!(I_idx, i)
				push!(J_idx, j)
				push!(V_val, val)
				if i != j
					push!(I_idx, j)
					push!(J_idx, i)
					push!(V_val, val)
				end
			end
		end
	end

	println("Constructing sparse matrix from $(length(V_val)) triplets...")
	@time CI_Matrix = sparse(I_idx, J_idx, V_val, Ndet, Ndet)
	println("Sparse matrix construction complete.")
	println("Diagonalizing CI Matrix ...")
	n_roots = min(Ndet, 30)
	@time E_CI_raw, C_CI = Davidson(CI_Matrix, n_roots)
	println("CI Matrix diagonalization complete.")

	Ee_CI = E_CI_raw[1]
	VNN_CI = SCF_Results.VNN
	Etot_CI = Ee_CI + VNN_CI
	Ecorr = Etot_CI - SCF_Results.Etot
	EStates = E_CI_raw .+ VNN_CI

	@printf("\n--- Final Energy Results (CI) ---\n")
	@printf("Total Energy      = %.10f Hartree\n", Etot_CI)
	@printf("Electronic Energy = %.10f Hartree\n", Ee_CI)
	@printf("Nuclear Repulsion =  %.10f Hartree\n", VNN_CI)

	Dets_Vec = [unpack_det(d) for d in Dets]

	return CIResults(
		SCF_Results.Etot,
		Etot_CI,
		Ecorr,
		EStates,
		C_CI,
		Dets_Vec,
		MaxExcitation,
		Ndet,
	)
end

function RunUCI(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int)
	TStart = time_ns()
	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")

	SCF_Results = UHF_SCF(Molecule, Charge, Multiplicity, MaxIter = 128, Threshold = 1e-8)
	if isnothing(SCF_Results)
		error("UHF calculation did not converge. Aborting.")
		return
	end

	CI_Results = CalcUCI(SCF_Results, MaxExcitation)
	if isnothing(CI_Results)
		error("CI calculation did not converge. Aborting.")
		return
	end

	println("\n--- Post-CI Analysis ---")
	println("\nGround State Wavefunction Analysis (Top 5 components):")
	GroundStateVector = CI_Results.CIVectors[:, 1]
	sorted_indices = sortperm(abs.(GroundStateVector), rev = true)
	RefDeterminant = CI_Results.DeterminantSpace[1]

	@printf("  Coeff    | Contribution  | Determinant Configuration\n")
	println("  -----------------------------------------------------")
	for i in 1:min(5, CI_Results.NumberOfDeterminants)
		idx = sorted_indices[i]
		coeff = GroundStateVector[idx]
		determinant = CI_Results.DeterminantSpace[idx]
		det_str = string(determinant)
		is_ref = (determinant == RefDeterminant) ? "(Ref)" : ""
		@printf("  %-8.4f | %-12.2f%% | %s %s\n", coeff, coeff^2 * 100, det_str, is_ref)
	end
	println("  -----------------------------------------------------\n")

	TEnd = time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")

	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia UHF-CI at $(DateTime).")
	return (Etot = CI_Results.EtotCI,)
end

function RunRCI(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int)
	TStart = time_ns()
	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")

	RHF_Res = RHF_SCF(Molecule, Charge, Multiplicity; MaxIter = 128, Threshold = 1e-8)
	if isnothing(RHF_Res)
		error("RHF calculation did not converge. Aborting.")
		return
	end
	SCF_Results = RHF2UHF(RHF_Res)

	CI_Results = CalcUCI(SCF_Results, MaxExcitation)
	if isnothing(CI_Results)
		error("CI calculation did not converge. Aborting.")
		return
	end

	println("\n--- Post-CI Analysis ---")
	println("\nGround State Wavefunction Analysis (Top 5 components):")
	GroundStateVector = CI_Results.CIVectors[:, 1]
	sorted_indices = sortperm(abs.(GroundStateVector), rev = true)
	RefDeterminant = CI_Results.DeterminantSpace[1]

	@printf("  Coeff    | Contribution  | Determinant Configuration\n")
	println("  -----------------------------------------------------")
	for i in 1:min(5, CI_Results.NumberOfDeterminants)
		idx = sorted_indices[i]
		coeff = GroundStateVector[idx]
		determinant = CI_Results.DeterminantSpace[idx]
		det_str = string(determinant)
		is_ref = (determinant == RefDeterminant) ? "(Ref)" : ""
		@printf("  %-8.4f | %-12.2f%% | %s %s\n", coeff, coeff^2 * 100, det_str, is_ref)
	end
	println("  -----------------------------------------------------\n")

	TEnd = time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")

	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia RHF-CI at $(DateTime).")
	return (Etot = CI_Results.EtotCI,)
end
