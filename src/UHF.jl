export RunUHF, UHFResults, UHF_SCF

struct UHFResults
	Molecule::Vector{Atom}
	BasisSet::Vector{Basis}
	ENumAlpha::Int
	ENumBeta::Int
	S::Matrix{Float64}
	Hcore::Matrix{Float64}
	ERI::Array{Float64, 4}
	CAlpha::Matrix{Float64}
	CBeta::Matrix{Float64}
	EAlpha::Vector{Float64}
	EBeta::Vector{Float64}
	Ee::Float64
	VNN::Float64
	Etot::Float64
end





function UHF_SCF(Molecule::Vector{Atom}, charge::Int, multiplicity::Int; MaxIter = 128, Threshold = 1e-7)
	BasisSet = generate_basis_list(Molecule)
	BNum = length(BasisSet)
	ENum = sum(atom.Z for atom in Molecule) - charge
	ENumAlpha = (ENum + multiplicity - 1) ÷ 2
	ENumBeta = (ENum - multiplicity + 1) ÷ 2

	println("--- System Information ---")
	@printf("Basis functions: %d\n", BNum)
	@printf("Alpha electrons: %d\n", ENumAlpha)
	@printf("Beta electrons:  %d\n", ENumBeta)
	println("--------------------------\n")

	S, T, V, ERI = CalcMatrices(BasisSet, Molecule)
	Hcore = T + V
	X = S^(-0.5)
	EGuess, CGuess = eigen(X' * Hcore * X)
	C = X * CGuess[:, sortperm(EGuess)]
	PAlpha = C[:, 1:ENumAlpha] * C[:, 1:ENumAlpha]'
	PBeta = C[:, 1:ENumBeta] * C[:, 1:ENumBeta]'

	DIIS = DIISManager{Tuple{Matrix{Float64}, Matrix{Float64}}}(10)

	println("\n--- Starting SCF Iterations ---")
	Etot_old = 0.0

	for i in 1:MaxIter
		Ptotal = PAlpha + PBeta
		GAlpha = [sum(Ptotal[k, l]*ERI[i, j, k, l] - PAlpha[k, l]*ERI[i, l, k, j] for k in 1:BNum, l in 1:BNum) for i in 1:BNum, j in 1:BNum]
		GBeta = [sum(Ptotal[k, l]*ERI[i, j, k, l] - PBeta[k, l]*ERI[i, l, k, j] for k in 1:BNum, l in 1:BNum) for i in 1:BNum, j in 1:BNum]
		FAlpha_current = Hcore + GAlpha
		FBeta_current = Hcore + GBeta

		ErrMat = X' * (FAlpha_current * PAlpha * S - S * PAlpha * FAlpha_current) * X +
				 X' * (FBeta_current * PBeta * S - S * PBeta * FBeta_current) * X

		ZeroMat = zeros(Float64, size(ErrMat))
		FAlpha, FBeta = diis_update!(DIIS, (FAlpha_current, FBeta_current), (ErrMat, ZeroMat))

		EAlphaVec, CprimeAlpha = eigen(X' * FAlpha * X)
		EBetaVec, CprimeBeta = eigen(X' * FBeta * X)

		pAlpha = sortperm(EAlphaVec)
		pBeta = sortperm(EBetaVec)
		EAlphaVec = EAlphaVec[pAlpha]
		EBetaVec = EBetaVec[pBeta]

		CAlpha = X * CprimeAlpha[:, sortperm(EAlphaVec)]
		CBeta = X * CprimeBeta[:, sortperm(EBetaVec)]

		PnewAlpha = CAlpha[:, 1:ENumAlpha] * CAlpha[:, 1:ENumAlpha]'
		PnewBeta = CBeta[:, 1:ENumBeta] * CBeta[:, 1:ENumBeta]'

		VNN = sum(Molecule[i].Z * Molecule[j].Z / norm(Molecule[i].position .- Molecule[j].position) for i in 1:length(Molecule) for j in (i+1):length(Molecule))
		Ee = 0.5 * sum(PAlpha .* (Hcore + FAlpha)) + 0.5 * sum(PBeta .* (Hcore + FBeta))
		Etot = Ee + VNN
		delta_E = abs(Etot - Etot_old)
		delta_P = max(sqrt(sum((PnewAlpha - PAlpha) .^ 2)), sqrt(sum((PnewBeta - PBeta) .^ 2)))
		@printf("Iteration %3d: E = %-16.10f  ΔE = %-12.2e  ΔP = %.2e\n", i, Etot, delta_E, delta_P)
		PAlpha = PnewAlpha
		PBeta = PnewBeta
		Etot_old = Etot
		if delta_E < Threshold && delta_P < Threshold
			println("\nSCF converged in $i iterations.")
			@printf("\n--- Final Energy Results ---\n")
			@printf("Electronic Energy = %.10f Hartree\n", Ee)
			@printf("Nuclear Repulsion = %.10f Hartree\n", VNN)
			@printf("Total Energy      = %.10f Hartree\n", Etot)
			println("----------------------------\n")
			return UHFResults(Molecule, BasisSet, ENumAlpha, ENumBeta, S, Hcore, ERI, CAlpha, CBeta, EAlphaVec, EBetaVec, Ee, VNN, Etot)
		end
	end
	println("\nSCF failed to converge after $MaxIter iterations.")
	return nothing
end


function RunUHF(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int; MaxIter = 100, Threshold = 1e-8)
	TStart=time_ns()
	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")
	SCF_Results=UHF_SCF(Molecule, Charge, Multiplicity, MaxIter = MaxIter, Threshold = Threshold)
	if isnothing(SCF_Results)
		error("UHF calculation did not converge. Aborting.")
		return
	end
	TEnd=time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia UHF at $(DateTime).")
	return SCF_Results
end
