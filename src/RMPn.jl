antisym(G_MO, p, q, r, s) = G_MO[p, q, r, s] - G_MO[p, q, s, r]

export RunRMPn

function TransG(G_AO::Array{Float64, 4}, C::Matrix{Float64})
    NBasis = size(C, 1) 
    NMO = size(C, 2)
    println("\n--- Starting Integral Transformation ---")
    G_tmp1 = zeros(NMO, NBasis * NBasis * NBasis)
    mul!(G_tmp1, C', reshape(G_AO, (NBasis, NBasis * NBasis * NBasis)))
    tmp_tensor_1 = reshape(G_tmp1, (NMO, NBasis, NBasis, NBasis))
    permuted_1 = permutedims(tmp_tensor_1, (2, 3, 4, 1))
    G_tmp2 = zeros(NMO, NBasis * NBasis * NMO)
    mul!(G_tmp2, C', reshape(permuted_1, (NBasis, NBasis * NBasis * NMO)))
    tmp_tensor_2 = reshape(G_tmp2, (NMO, NBasis, NBasis, NMO))
    permuted_2 = permutedims(tmp_tensor_2, (2, 3, 4, 1))
    G_tmp3 = zeros(NMO, NBasis * NMO * NMO)
    mul!(G_tmp3, C', reshape(permuted_2, (NBasis, NBasis * NMO * NMO)))
    tmp_tensor_3 = reshape(G_tmp3, (NMO, NBasis, NMO, NMO))
    permuted_3 = permutedims(tmp_tensor_3, (2, 3, 4, 1))
    G_MO_flat = zeros(NMO, NMO * NMO * NMO)
    mul!(G_MO_flat, C', reshape(permuted_3, (NBasis, NMO * NMO * NMO)))
    G_MO_perm = reshape(G_MO_flat, (NMO, NMO, NMO, NMO))
    G_MO = permutedims(G_MO_perm, (2, 3, 4, 1))
    println("--- Integral Transformation Completed ---\n")
    return G_MO
end

function CalcRMP2(G_MO::Array{Float64, 4}, E_MO::Vector{Float64}, NOcc::Int)
	NMO = length(E_MO)
	NVir = NMO - NOcc

	EOcc = E_MO[1:NOcc]
	EVir = E_MO[(NOcc+1):NMO]
	TAmp = zeros(NOcc, NOcc, NVir, NVir)

	E_RMP2_Corr = 0.0

	for i in 1:NOcc, j in 1:NOcc
		for a in 1:NVir, b in 1:NVir
			aIdx = a + NOcc
			bIdx = b + NOcc

			Denom = EOcc[i] + EOcc[j] - EVir[a] - EVir[b]

			integral_iajb = G_MO[i, aIdx, j, bIdx]
			integral_ibja = G_MO[i, bIdx, j, aIdx]

			TAmp[i, j, a, b] = (2 * integral_iajb - integral_ibja) / Denom

			E_RMP2_Corr += integral_iajb * TAmp[i, j, a, b]
		end
	end

	return E_RMP2_Corr, TAmp
end

function RMPn(RHF_Results::RHFResults, order::Int)
	println("--- Starting Restricted MPn Calculation ---")

	C = RHF_Results.C
	e = RHF_Results.E
	G_AO = RHF_Results.ERI
	ENum = RHF_Results.ENum
	EtotRHF = RHF_Results.Etot

	NOcc = ENum ÷ 2
	G_MO = TransG(G_AO, C)

	E_RMP2_Corr, TAmp = CalcRMP2(G_MO, e, NOcc)
	E_RMP2_total = EtotRHF + E_RMP2_Corr

	@printf("\n--- MP2 Energy Results ---\n")
	@printf("MP2 Correlation Energy = %.10f Hartree\n", E_RMP2_Corr)
	@printf("Total MP2 Energy       = %.10f Hartree\n", E_RMP2_total)
	println("----------------------------\n")
	if order==2
		return (E_RMP2_Corr = E_RMP2_Corr, E_RMP2_total = E_RMP2_total)
	else
		println("Higher order MPn not implemented yet.")
		return nothing
	end
end

function RunRMPn(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, order::Int)
	TStart=time_ns()

	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	RHF_Results = RHF_SCF(Molecule, Charge, Multiplicity, MaxIter = 100, Threshold = 1e-8)

	if RHF_Results !== nothing
		println("RHF Calculation Successful. Proceeding to MPn.")

		MPn_Results = RMPn(RHF_Results, order)
		println("\n--- Calculation Summary ---")
		println("RHF Total Energy: ", RHF_Results.Etot, " Hartree")
		if MPn_Results !== nothing
			println("MP2 Total Energy: ", MPn_Results.E_RMP2_total, " Hartree")
			TEnd=time_ns()
			TSeconds = (TEnd - TStart) / 1e9
			days = floor(Int, TSeconds / 86400)
			hours = floor(Int, (TSeconds % 86400) / 3600)
			minutes = floor(Int, (TSeconds % 3600) / 60)
			seconds = TSeconds % 60
			DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
			@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
			@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
			println(" Normal termination of Julia RMPn at $(DateTime).")
		end
		println("-------------------------\n")
		return (Etot = MPn_Results.E_RMP2_total,)
	else
		println("RHF calculation failed to converge. Cannot perform MPn calculation.")
	end
end