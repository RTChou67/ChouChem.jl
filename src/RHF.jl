using ChemAlgebra: DIISManager, diis_update!

struct RHFResults
	Molecule::Vector{Atom}
	BasisSet::Vector{Basis}
	ENum::Int
	S::Matrix{Float64}
	Hcore::Matrix{Float64}
	ERI::Array{Float64, 4}
	C::Matrix{Float64}
	E::Vector{Float64}
	P::Matrix{Float64}
	Ee::Float64
	VNN::Float64
	Etot::Float64
end

function RHF_SCF(Molecule::Vector{Atom}, Charge::Int,	Multiplicity::Int; MaxIter = 100, Threshold = 1e-8)
	BasisSet = generate_basis_list(Molecule)
	BNum = length(BasisSet)
	ENum = sum(atom.Z for atom in Molecule) - Charge
	Nocc = ENum ÷ 2
	println("--- System Information ---")
	@printf("Basis functions: %d\n", BNum)
	@printf("Electrons:       %d\n", ENum)
	println("--------------------------\n")
	S, T, V, ERI = CalcMatrices(BasisSet, Molecule)
	Hcore = T + V
	X = S^(-0.5)
	E_guess, C_guess = eigen(X' * Hcore * X)
	p = sortperm(E_guess)
	C = X * C_guess[:, p]
	P = 2 * C[:, 1:Nocc] * C[:, 1:Nocc]'
	DIIS = DIISManager{Matrix{Float64}}(10)
	println("--- Starting SCF Iterations (with DIIS) ---")
	Etot_old = 0.0
	for i in 1:MaxIter
		G = [sum(P[k, l] * (ERI[i, j, k, l] - 0.5 * ERI[i, l, k, j]) for k in 1:BNum, l in 1:BNum) for i in 1:BNum, j in 1:BNum]
		F_current = Hcore + G
		ErrMat = X' * (F_current * P * S - S * P * F_current) * X
		
		F = diis_update!(DIIS, F_current, ErrMat)
		
		Fprime = X' * F * X
		E, Cprime = eigen(Fprime)
		p = sortperm(E)
		E = E[p]
		C = X * Cprime[:, p]
		Pnew = 2 * C[:, 1:Nocc] * C[:, 1:Nocc]'
		VNN = sum(Molecule[i].Z * Molecule[j].Z / norm(Molecule[i].position .- Molecule[j].position) for i in 1:length(Molecule) for j in (i+1):length(Molecule))
		Ee = 0.5 * sum(P .* (Hcore + F))
		Etot = Ee + VNN
		delta_E = abs(Etot - Etot_old)
		delta_P = sqrt(sum((Pnew - P) .^ 2))
		@printf("Iteration %3d: E = %-16.10f  ΔE = %-12.2e  ΔP = %.2e\n", i, Etot, delta_E, delta_P)
		P = Pnew
		Etot_old = Etot
		if delta_E < Threshold && delta_P < Threshold
			println("\nSCF converged in $i iterations.")
			@printf("\n--- Final Energy Results ---\n")
			@printf("Electronic Energy = %.10f Hartree\n", Ee)
			@printf("Nuclear Repulsion = %.10f Hartree\n", VNN)
			@printf("Total Energy      = %.10f Hartree\n", Etot)
			println("----------------------------\n")
			return RHFResults(Molecule, BasisSet, ENum, S, Hcore, ERI, C, E, P, Ee, VNN, Etot)
		end
	end
	println("\nSCF failed to converge after $MaxIter iterations.")
	return nothing
end

function RunRHF(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int)
	TStart = time_ns()
	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")
	SCF_Results = RHF_SCF(Molecule, Charge, Multiplicity, MaxIter = 100, Threshold = 1e-8)
	if isnothing(SCF_Results)
		error("RHF calculation did not converge. Aborting.")
		return
	end
	TEnd = time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia RHF at $(DateTime).")
	return SCF_Results
end