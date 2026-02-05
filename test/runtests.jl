using Test
using ChouChem
using Printf

function mkMol(Coord, BasisSet::String)
	Mol = Vector{Atom}()
	for atom in Coord
		push!(Mol, Atom(atom[1], atom[2], BasisSet, atom[3]))
	end
	return Mol
end

HFCoord = [
	("H", 1, (0.000, 0.000, 0.000)),
	("F", 9, (0.000, 0.000, 1.500)),
]

H2OCoord = [
	("O", 8, (0.0000, 0.0000, 0.0000)),
	("H", 1, (0.0000, 1.0000, 0.8000)),
	("H", 1, (0.0000, -1.0000, 0.8000)),
]

HFSTO3G = mkMol(HFCoord, "STO-3G")
HF631G = mkMol(HFCoord, "6-31G")
H2OSTO3G = mkMol(H2OCoord, "STO-3G")
H2O631G = mkMol(H2OCoord, "6-31G")

# Structure to hold results
TestResults = Vector{Tuple{String, String, String, Float64}}()

@testset "ChouChem Tests" begin
	@testset "HF Molecule" begin
		@testset "STO-3G" begin
			@testset "RHF" begin
				HF_RHF_Results = RunRHF(HFSTO3G, 0, 1)
				@test HF_RHF_Results !== nothing
				push!(TestResults, ("RHF", "HF", "STO-3G", HF_RHF_Results.Etot))
			end
			@testset "UHF" begin
				HF_UHF_Results = RunUHF(HFSTO3G, 1, 2; MaxIter=512)
				@test HF_UHF_Results !== nothing
				push!(TestResults, ("UHF", "HF+", "STO-3G", HF_UHF_Results.Etot))
			end
			@testset "RMP2" begin
				HF_RMP2_Results = RunRMPn(HFSTO3G, 0, 1, 2)
				@test HF_RMP2_Results !== nothing
				push!(TestResults, ("RMP2", "HF", "STO-3G", HF_RMP2_Results.Etot))
			end
			@testset "UCI" begin
				HF_UCI_Results = RunUCI(HFSTO3G, 1, 2, 2; MaxIter=512)
				@test HF_UCI_Results !== nothing
				push!(TestResults, ("UCI", "HF+", "STO-3G", HF_UCI_Results.Etot))
			end
			@testset "RCI" begin
				HF_RCI_Results = RunRCI(HFSTO3G, 0, 1, 2)
				@test HF_RCI_Results !== nothing
				push!(TestResults, ("RCI", "HF", "STO-3G", HF_RCI_Results.Etot))
			end
		end
		@testset "6-31G" begin
			@testset "RHF" begin
				HF_RHF_Results = RunRHF(HF631G, 0, 1)
				@test HF_RHF_Results !== nothing
				push!(TestResults, ("RHF", "HF", "6-31G", HF_RHF_Results.Etot))
			end
			@testset "UHF" begin
				HF_UHF_Results = RunUHF(HF631G, 1, 2; MaxIter=512)
				@test HF_UHF_Results !== nothing
				push!(TestResults, ("UHF", "HF+", "6-31G", HF_UHF_Results.Etot))
			end
			@testset "RMP2" begin
				HF_RMP2_Results = RunRMPn(HF631G, 0, 1, 2)
				@test HF_RMP2_Results !== nothing
				push!(TestResults, ("RMP2", "HF", "6-31G", HF_RMP2_Results.Etot))
			end
			@testset "UCI" begin
				HF_UCI_Results = RunUCI(HF631G, 1, 2, 2; MaxIter=512)
				@test HF_UCI_Results !== nothing
				push!(TestResults, ("UCI", "HF+", "6-31G", HF_UCI_Results.Etot))
			end
			@testset "RCI" begin
				HF_RCI_Results = RunRCI(HF631G, 0, 1, 2)
				@test HF_RCI_Results !== nothing
				push!(TestResults, ("RCI", "HF", "6-31G", HF_RCI_Results.Etot))
			end
		end
	end
	@testset "H2O Molecule" begin
		@testset "STO-3G" begin
			@testset "RHF" begin
				H2O_RHF_Results = RunRHF(H2OSTO3G, 0, 1)
				@test H2O_RHF_Results !== nothing
				push!(TestResults, ("RHF", "H2O", "STO-3G", H2O_RHF_Results.Etot))
			end
			@testset "UHF" begin
				H2O_UHF_Results = RunUHF(H2OSTO3G, 1, 2; MaxIter=512)
				@test H2O_UHF_Results !== nothing
				push!(TestResults, ("UHF", "H2O+", "STO-3G", H2O_UHF_Results.Etot))
			end
			@testset "RMP2" begin
				H2O_RMP2_Results = RunRMPn(H2OSTO3G, 0, 1, 2)
				@test H2O_RMP2_Results !== nothing
				push!(TestResults, ("RMP2", "H2O", "STO-3G", H2O_RMP2_Results.Etot))
			end
			@testset "UCI" begin
				H2O_UCI_Results = RunUCI(H2OSTO3G, 1, 2, 2; MaxIter=512)
				@test H2O_UCI_Results !== nothing
				push!(TestResults, ("UCI", "H2O+", "STO-3G", H2O_UCI_Results.Etot))
			end
			@testset "RCI" begin
				H2O_RCI_Results = RunRCI(H2OSTO3G, 0, 1, 2)
				@test H2O_RCI_Results !== nothing
				push!(TestResults, ("RCI", "H2O", "STO-3G", H2O_RCI_Results.Etot))
			end
		end
		@testset "6-31G" begin
			@testset "RHF" begin
				H2O_RHF_Results = RunRHF(H2O631G, 0, 1)
				@test H2O_RHF_Results !== nothing
				push!(TestResults, ("RHF", "H2O", "6-31G", H2O_RHF_Results.Etot))
			end
			@testset "UHF" begin
				H2O_UHF_Results = RunUHF(H2O631G, 1, 2; MaxIter=512)
				@test H2O_UHF_Results !== nothing
				push!(TestResults, ("UHF", "H2O+", "6-31G", H2O_UHF_Results.Etot))
			end
			@testset "RMP2" begin
				H2O_RMP2_Results = RunRMPn(H2O631G, 0, 1, 2)
				@test H2O_RMP2_Results !== nothing
				push!(TestResults, ("RMP2", "H2O", "6-31G", H2O_RMP2_Results.Etot))
			end
			@testset "UCI" begin
				H2O_UCI_Results = RunUCI(H2O631G, 1, 2, 2; MaxIter=512)
				@test H2O_UCI_Results !== nothing
				push!(TestResults, ("UCI", "H2O+", "6-31G", H2O_UCI_Results.Etot))
			end
			@testset "RCI" begin
				H2O_RCI_Results = RunRCI(H2O631G, 0, 1, 2)
				@test H2O_RCI_Results !== nothing
				push!(TestResults, ("RCI", "H2O", "6-31G", H2O_RCI_Results.Etot))
			end
		end
	end

	@testset "RCCSD H-F Experiment" begin
		# Experimental bond length of H-F is ~0.9168 A
		HF_Exp = mkMol([
			("H", 1, (0.0, 0.0, 0.0)),
			("F", 9, (0.0, 0.0, 0.9168))
		], "STO-3G") # Using STO-3G for speed in test, user asked for functionality
		
		CCSD_Res = RunRCCSD(HF_Exp, 0, 1)
		@test CCSD_Res !== nothing
		@printf("RCCSD Correlation Energy: %.10f\n", CCSD_Res.Ecorr)
		push!(TestResults, ("RCCSD", "HF(Exp)", "STO-3G", CCSD_Res.EtotCCSD))
	end

	@testset "UCCSD H-F Experiment" begin
		HF_Exp = mkMol([
			("H", 1, (0.0, 0.0, 0.0)),
			("F", 9, (0.0, 0.0, 0.9168))
		], "STO-3G")
		
		UCCSD_Res = RunUCCSD(HF_Exp, 1, 2; MaxIter=512)
		@test UCCSD_Res !== nothing
		@printf("UCCSD Correlation Energy: %.10f\n", UCCSD_Res.Ecorr)
		push!(TestResults, ("UCCSD", "HF+(Exp)", "STO-3G", UCCSD_Res.EtotCCSD))
	end
end

println("\n" * "="^70)
println(" "^25 * "TEST RESULTS SUMMARY")
println("="^70)
@printf("%-10s | %-10s | %-10s | %s\n", "Method", "Molecule", "Basis", "Total Energy (Hartree)")
println("-"^70)
for r in TestResults
	@printf("%-10s | %-10s | %-10s | %.10f\n", r[1], r[2], r[3], r[4])
end
println("="^70 * "\n")
