export CalcMatrices, RHF2UHF

function CalcMatrices(BasisSet, Molecule)
	BNum = length(BasisSet)
	TimeS1=time_ns()
	S = [Sij(BasisSet[i], BasisSet[j]) for i in 1:BNum, j in 1:BNum]
	TimeS2=time_ns()
	println("Calculation for S Matrix took $( (TimeS2-TimeS1)/1e6 ) ms")
	TimeT1=time_ns()
	T = [Tij(BasisSet[i], BasisSet[j]) for i in 1:BNum, j in 1:BNum]
	TimeT2=time_ns()
	println("Calculation for T Matrix took $( (TimeT2-TimeT1)/1e6 ) ms")
	TimeV1=time_ns()
	V = [Vij(BasisSet[i], BasisSet[j], Molecule) for i in 1:BNum, j in 1:BNum]
	TimeV2=time_ns()
	println("Calculation for V Matrix took $( (TimeV2-TimeV1)/1e6 ) ms")
	TimeERI1=time_ns()
	ERI = [Gijkl(BasisSet[i], BasisSet[j], BasisSet[k], BasisSet[l]) for i in 1:BNum, j in 1:BNum, k in 1:BNum, l in 1:BNum]
	TimeERI2=time_ns()
	println("Calculation for ERI Tensor took $( (TimeERI2-TimeERI1)/1e6 ) ms")
	return S, T, V, ERI
end

function RHF2UHF(rhf::RHFResults)
	return UHFResults(
		rhf.Molecule,
		rhf.BasisSet,
		rhf.ENum ÷ 2,
		rhf.ENum ÷ 2,
		rhf.S,
		rhf.Hcore,
		rhf.ERI,
		rhf.C,
		rhf.C,
		rhf.E,
		rhf.E,
		rhf.Ee,
		rhf.VNN,
		rhf.Etot,
	)
end

function pack_det(indices::Vector{Int})::Int
	d = zero(Int)
	for i in indices
		d |= (one(Int) << (i - 1))
	end
	return d
end

function unpack_det(d::Int)::Vector{Int}
	res = Int[]
	idx = 1
	while d > 0
		if (d & 1) == 1
			push!(res, idx)
		end
		d >>= 1
		idx += 1
	end
	return res
end

@inline function get_position(det::Int, k::Int)
	mask = (one(Int) << (k - 1)) - 1
	return count_ones(det & mask)
end

function get_phase_single(det_i::Int, det_j::Int, p::Int, q::Int)
	pos_p = get_position(det_i, p)
	pos_q = get_position(det_j, q)
	return iseven(pos_p + pos_q) ? 1.0 : -1.0
end

function get_phase_double(det_i::Int, det_j::Int, p::Int, q::Int, r::Int, s::Int)
	pos_p = get_position(det_i, p)
	pos_q = get_position(det_i, q)

	pos_r = get_position(det_j, r)
	pos_s = get_position(det_j, s)

	return iseven(pos_p + pos_q + pos_r + pos_s) ? 1.0 : -1.0
end

function TransERI_blas(ERI_AO::Array{Float64, 4}, c1::Matrix{Float64}, c2::Matrix{Float64}, c3::Matrix{Float64}, c4::Matrix{Float64}, ONum::Int)
	N = ONum
	M_AO = reshape(ERI_AO, (N, N*N*N))
	M_tmp1 = zeros(N, N*N*N)
	mul!(M_tmp1, c1', M_AO)

	tmp_pjkl = reshape(M_tmp1, (N, N, N, N))
	tmp_pjkl_perm = permutedims(tmp_pjkl, (1, 3, 4, 2))
	M_tmp1_perm = reshape(tmp_pjkl_perm, (N*N*N, N))
	M_tmp2 = zeros(N*N*N, N)
	mul!(M_tmp2, M_tmp1_perm, c2)

	tmp_pqkl_perm = reshape(M_tmp2, (N, N, N, N))
	tmp_pqkl = permutedims(tmp_pqkl_perm, (1, 4, 2, 3))
	tmp_pqkl_perm = permutedims(tmp_pqkl, (1, 2, 4, 3))
	M_tmp2_perm = reshape(tmp_pqkl_perm, (N*N*N, N))
	M_tmp3 = zeros(N*N*N, N)
	mul!(M_tmp3, M_tmp2_perm, c3)

	tmp_pqrl_perm = reshape(M_tmp3, (N, N, N, N))
	tmp_pqrl = permutedims(tmp_pqrl_perm, (1, 2, 4, 3))
	M_tmp3_perm = reshape(tmp_pqrl, (N*N*N, N))
	ERI_MO_M = zeros(N*N*N, N)
	mul!(ERI_MO_M, M_tmp3_perm, c4)
	return reshape(ERI_MO_M, (N, N, N, N))
end

function build_spin_eri(aaaa, bbbb, aabb, ONum::Int)
	N = ONum
	N2 = 2 * N
	ERI = zeros(Float64, N2, N2, N2, N2)
	ERI[1:N, 1:N, 1:N, 1:N] .= aaaa
	ERI[(N+1):2N, (N+1):2N, (N+1):2N, (N+1):2N] .= bbbb
	ERI[1:N, 1:N, (N+1):2N, (N+1):2N] .= aabb
	ERI[(N+1):2N, (N+1):2N, 1:N, 1:N] .= permutedims(aabb, (3, 4, 1, 2))
	return ERI
end