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
Molecules = [
	("HF", HFCoord),
	("H2O", H2OCoord)
]

BasisSets = ["STO-3G", "6-31G"]
TestConfigs = [
	(0, 1, [
		("RHF",   RunRHF,   (),                 "RHF"),
		("RMP2",  RunRMPn,  (2,),               "MP2"),
		("RCI",   RunRCI,   (2,),               "CISD"),
		("RCCSD", RunRCCSD, (),                 "CCSD")
	]),
	(1, 2, [
		("UHF",   RunUHF,   (),                 "UHF"),
		("UCI",   RunUCI,   (2,),               "CISD"),
		("UCCSD", RunUCCSD, (),                 "CCSD")
	])
]

TestResults = Vector{Tuple{String, String, String, String, Float64}}() # Mol, Basis, Charge, Method, Energy
G16_Dir = "test_g16"

GaussianRef = Dict{Tuple{String, String, String, String}, Float64}(
    ("HF", "STO-3G", "0", "RHF") => -98.4147505,
    ("HF", "STO-3G", "0", "RMP2") => -98.4731522,
    ("HF", "STO-3G", "0", "RCI") => -98.5193602,
    ("HF", "STO-3G", "0", "RCCSD") => -98.5193602,
    ("HF", "STO-3G", "1", "UHF") => -98.1167715,
    ("HF", "STO-3G", "1", "UCI") => -98.1334809,
    ("HF", "STO-3G", "1", "UCCSD") => -98.1416672,
    ("HF", "6-31G", "0", "RHF") => -99.8550550,
    ("HF", "6-31G", "0", "RMP2") => -100.0142173,
    ("HF", "6-31G", "0", "RCI") => -100.0142910,
    ("HF", "6-31G", "0", "RCCSD") => -100.0226417,
    ("HF", "6-31G", "1", "UHF") => -99.4226454,
    ("HF", "6-31G", "1", "UCI") => -99.5203610,
    ("HF", "6-31G", "1", "UCCSD") => -99.5226893,
    ("H2O", "STO-3G", "0", "RHF") => -74.8492756,
    ("H2O", "STO-3G", "0", "RMP2") => -74.9219032,
    ("H2O", "STO-3G", "0", "RCI") => -74.9522844,
    ("H2O", "STO-3G", "0", "RCCSD") => -74.9568665,
    ("H2O", "STO-3G", "1", "UHF") => -74.5963826,
    ("H2O", "STO-3G", "1", "UCI") => -74.6657236,
    ("H2O", "STO-3G", "1", "UCCSD") => -74.6676946,
    ("H2O", "6-31G", "0", "RHF") => -75.8727853,
    ("H2O", "6-31G", "0", "RMP2") => -76.0331925,
    ("H2O", "6-31G", "0", "RCI") => -76.0327560,
    ("H2O", "6-31G", "0", "RCCSD") => -76.0427793,
    ("H2O", "6-31G", "1", "UHF") => -75.5120875,
    ("H2O", "6-31G", "1", "UCI") => -75.633082,
    ("H2O", "6-31G", "1", "UCCSD") => -75.6384025
)

@testset "ChouChem Tests" begin
	for (mol_name, coords) in Molecules
		@testset "$mol_name" begin
			for basis in BasisSets
				@testset "$basis" begin
					MolObj = mkMol(coords, basis)

					for (charge, mult, methods) in TestConfigs
						charge_str = charge == 0 ? "Neutral" : "Cation"
						@testset "$charge_str (Q=$charge, M=$mult)" begin
							for (method_name, func, args, g16_key) in methods
								@testset "$method_name" begin
									
									local result = nothing
									try
										if method_name == "RHF"
											result = func(MolObj, charge, mult)
										elseif method_name == "UHF"
											result = func(MolObj, charge, mult; MaxIter=512)
										elseif method_name == "RMP2"
											result = func(MolObj, charge, mult, args...)
										elseif method_name == "UCI"
											result = func(MolObj, charge, mult, args...; MaxIter=512)
										elseif method_name == "RCI"
											result = func(MolObj, charge, mult, args...)
										elseif method_name == "RCCSD"
											result = func(MolObj, charge, mult)
										elseif method_name == "UCCSD"
											result = func(MolObj, charge, mult; MaxIter=512)
									end
								catch e
									println("Error running $method_name for $mol_name/$basis: $e")
									throw(e)
								end

									@test result !== nothing
									
									energy = 0.0
									if hasproperty(result, :Etot)
										energy = result.Etot
									elseif hasproperty(result, :EtotMPn)
										energy = result.EtotMPn
									elseif hasproperty(result, :EtotCI)
										energy = result.EtotCI
									elseif hasproperty(result, :EtotCCSD)
										energy = result.EtotCCSD
									end
									
									push!(TestResults, (mol_name, basis, string(charge), method_name, energy))
								end
							end
						end
					end
				end
			end
		end
	end
end


println("\n" * "="^120)
println(" "^45 * "TEST RESULTS SUMMARY")
println("="^120)
@printf("%-8s | %-8s | %-3s | %-8s | %-16s | %-16s | %s\n", "Mol", "Basis", "Q", "Method", "ChouChem", "Gaussian", "Diff (Ha)")
println("-"^120)
for r in TestResults
    mol, basis, q, method, chou_e = r
    g16_e = get(GaussianRef, (mol, basis, q, method), NaN)
    diff = abs(chou_e - g16_e)
    status = isnan(g16_e) ? "N/A" : (diff < 1e-6 ? "OK" : "DIFF")
    if status == "DIFF" && diff < 1e-3
        status = "DIFF(Small)"
    elseif status == "DIFF" && diff < 5e-2
        status = "DIFF(Med)"
    end
    
	@printf("%-8s | %-8s | %-3s | %-8s | %16.10f | %16.10f | %10.2e %s\n", mol, basis, q, method, chou_e, g16_e, diff, status)
end
println("="^120 * "\n")
