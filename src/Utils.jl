export BOHR2ANG, ang2bohr, print_molecule, print_job_time, print_ci_analysis

const BOHR2ANG = 0.52917721092

function ang2bohr(MolInAng::Vector{Atom})
	return [Atom(a.symbol, a.Z, a.basis_set, a.position ./ BOHR2ANG) for a in MolInAng]
end

function print_molecule(MolInAng::Vector{Atom})
	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")
end

function print_job_time(TStart::UInt64, label::String)
	TSeconds = (time_ns() - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia $label at $(DateTime).")
end

function print_ci_analysis(CI_Results)
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
end
