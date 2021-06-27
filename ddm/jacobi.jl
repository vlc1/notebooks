### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 110d3d54-9225-4a44-83ce-3a7f509232fd
using OffsetArrays

# ╔═╡ 06626098-d059-11eb-0434-6f829bce295f
md"""
Let's first examine the Jacobi preconditioner with no overlap.

# Uniform partition of a rectangular array

```math
n_i = \left \lfloor \frac{n i}{p} \right \rfloor - \left \lfloor \frac{n \left ( i - 1 \right )}{p} \right \rfloor, \quad i \in \left \{ 1, \ldots, p \right \}.
```

"""

# ╔═╡ e6a47c50-0d83-4333-af8e-6b493aef773c
begin
	struct UniformPartition{N,T<:CartesianIndices{N}} <: AbstractArray{T,N}
		indices::T
		nover::NTuple{N,Int}
		nproc::NTuple{N,Int}
	end
end

# ╔═╡ 74880f2a-cd4b-49af-a574-c680785b9be6
getindices(p::UniformPartition) = p.indices

# ╔═╡ 8683a8d9-af47-4c3b-a075-b7b4796c297f
getover(p::UniformPartition) = p.nover

# ╔═╡ 5429496f-9473-48d2-9842-50bbace1403c
getproc(p::UniformPartition) = p.nproc

# ╔═╡ e630dee3-1781-4d81-9185-45096a87047c
Base.size(p::UniformPartition) = getproc(p)

# ╔═╡ 2b2cf2c0-c38d-4991-8a33-7a3b481754ec
Base.getindex(part::UniformPartition, index::CartesianIndex) =
	getindex(part, Tuple(index)...)

# ╔═╡ 087bc9b7-817d-49fc-ab00-513ab03f342a
function Base.getindex(part::UniformPartition, index...)
	indices = getindices(part)
	nover = getover(part)
	nproc = getproc(part)

	lo = Tuple(first(indices))
	hi = Tuple(last(indices))

	map(lo, hi, nover, nproc, index) do l, u, o, p, i
		d = u - l + 1
		UnitRange(
			max(d * (i - 1) ÷ p + l - o ÷ 2, l),
			min(d * i ÷ p - 1 + l + (o + 1) ÷ 2, u)
		)
	end |> CartesianIndices
end

# ╔═╡ f30e6e06-ca06-4a0e-8dc8-a9b45eb6db5c
partition = UniformPartition(
	# indices
	CartesianIndices((33:64, 17:32)),
	# overlap (nover)
	(2, 1),
	# processes (nproc)
	(9, 3)
)

# ╔═╡ a8d9eeff-bd77-4e02-89b0-3f676fa11fc3
md"""
# `CartesianIndices`

"""

# ╔═╡ 49bf2458-e192-4c01-80e4-f48b37bca3d2
function pad(indices, n)
	lo = Tuple(first(indices))
	hi = Tuple(last(indices))

	map(n, lo, hi) do i, l, h
		l - i:h + i
	end |> CartesianIndices
end

# ╔═╡ cb9a7522-8dba-4c95-a7c3-93ee8f2fd21b
md"""
# Ghosting

"""

# ╔═╡ baaa29dd-1fb6-451b-961c-20c9a45be294
begin
	struct GhostedArray{T,N,AA<:AbstractArray{T,N}} <: AbstractArray{T,N}
		parent::AA
		nghost::NTuple{N,Int}
	end

	function GhostedArray{T}(indices::CartesianIndices{N}, n) where {T,N}
		ghosted = pad(indices, n)
		data = OffsetArray(Array{T}(undef, size(ghosted)), ghosted)
		GhostedArray(data, n)
	end
end

# ╔═╡ d8452319-d87b-4d11-b21f-f2c8306179ee
Base.parent(a::GhostedArray) = a.parent

# ╔═╡ c72529da-9395-4409-9f55-b15e46611e0d
Base.getindex(a::GhostedArray, i...) = getindex(parent(a), i...)

# ╔═╡ a1a7383b-9585-4952-9a27-a6ed74b60ce0
Base.setindex!(a::GhostedArray, i...) = setindex!(parent(a), i...)

# ╔═╡ 3ecd7b07-e603-4994-8941-b704bae05299
# need to skip
Base.axes(a::GhostedArray) = axes(parent(a))

# ╔═╡ 0d2d52a4-c656-41dd-b641-17c209c3affa
Base.size(a::GhostedArray) = size(parent(a))

# ╔═╡ 2d8bf534-0203-4d58-9f2d-1aabd99176a3
nghost = (1, 1)

# ╔═╡ 0fecb93e-8d74-4c9d-a5b6-7fe43e175268
b = GhostedArray{Float64}(getindices(partition), nghost)

# ╔═╡ 449cd7fb-77f7-4c53-b665-0d322aed465d
axes(b)

# ╔═╡ d8298fed-3266-47c7-8b0a-1361a940c557
function initialize(::Type{T}, indices::CartesianIndices{N}, n) where {T,N}
	ghosted = pad(indices, n)
	OffsetArray(Array{T}(undef, size(ghosted)), ghosted)
end

# ╔═╡ a615ce0f-f57e-4e50-9334-4b487f9a2443
#b = initialize(Float64, getindices(partition), nst)

# ╔═╡ 37eb5392-9c17-4d06-883b-0d3693a954fc
view(b, getindices(partition))

# ╔═╡ b39d2736-c827-452b-a5a3-adbf19eb781d
Array{Float64}(undef, (2:3, 3:4))

# ╔═╡ 3c59a707-3b63-450d-90e1-1ed6c00da6c2
foo = begin
	local indices = ghostedoverlap(partition[1, 1], param...)
	OffsetArray(zeros(size(indices)), indices)
end

# ╔═╡ 86247385-1aed-4ed6-b70f-4a567201889e
bar = begin
	indices = semioverlap(partition[1, 1], param...)
	OffsetArray(view(foo, indices), indices)
end

# ╔═╡ cd8534b0-e522-4d42-8cc0-bf3ae75ad92e
#=
begin
	struct DDMArray{T,N,AA<:AbstractArray{T,N}} <: AbstractArray{T,N}
		parent::AA
    	nover::NTuple{N,Int}
		nst::NTuple{N,Int}

		# check parameters
	end

	DDMVector{T} = DDMArray{T,1}
end
=#

# ╔═╡ 47ce9c5c-4a82-43dc-b529-d83e55532e4c
#=
function Base.view(::Type{T}, a::DDMVector) where {T<:AbstractDomain}

	r = axes(T, parent(a), getoverlap(a), getstencil(a))

	s = first(r) + st + 2over:last(r) - st - 2over

	view(par, s)
end
=#

# ╔═╡ 0eba4d39-e3f7-45a1-96d6-847109e48e23
#overlap(a::DDMArray) = a.nover

# ╔═╡ 1f31f525-afc4-4e20-ae3e-fd3fb93aca24
#stencil(a::DDMArray) = a.nst

# ╔═╡ 46277571-c698-4dd7-a6a7-a6c609d0845e
#=
function init(nover, nst, interior)
	n = 2overlap + stencil
	[DDMArray(
		OffsetArray(
			Vector{Float64}(undef, length(r) + 2n),
			first(r) - n:last(r) + n),
		overlap,
		stencil) for r in interior]
end
=#

# ╔═╡ 5dab8c53-4fac-4ac4-833c-78918dd3d823
function initialize(x, partition, n)
	map(partition) do indices
		ghost = pad(indices, n)
		OffsetArray(similar(x, axes(ghost)), ghost)
	end
end

# ╔═╡ 5da69b07-3d42-4ee1-86ad-debee168de8d
#=
b = begin
	local indices = CartesianIndices(getdof(partition))
	local ghost = pad(indices, nst)
	OffsetArray(zeros(axes(ghost)), ghost)
end
=#

# ╔═╡ 1542f9dd-21e7-4e4e-a50e-5a234b1f8abe
CartesianIndices(getdof(partition))

# ╔═╡ dc293c44-578f-4d45-aca0-a6f766818314
bs = initialize(b, partition, nst)

# ╔═╡ f5d9153c-e31d-4562-a05b-19b4b7f0ca01
function update!(ys, x)
	for y in ys
		indices = CartesianIndices(y)
		y .= x[indices]
	end
end

# ╔═╡ 5182c5ca-b1f0-4773-97fe-33e9a939c540
bs[end]

# ╔═╡ 49118046-3d81-4e2c-80cf-0e438e3e391a
#bs = init(overlap, stencil, interior)

# ╔═╡ a00e640f-cec0-4a17-8215-a5d59e42fe07
#view(Interior, bs[1])

# ╔═╡ 798003a3-a0bc-47b8-b4ad-77a112da466e
#=
function update!(::Type{Interior}, ys, x)
	for y in ys
		interior = view(Interior, y)
		#interior .= x[]
	end
end
=#

# ╔═╡ 675b43e1-66ac-4143-8b9a-ae2ff3bad1d0
function update!(::Type{Ghost}, xs, interior)

	n = size(xs, 1)

	for (i, (center, r)) in enumerate(zip(xs, interior))

		overlap = getoverlap(center)
		stencil = getstencil(center)

		# left
		j = (n + i - 2) % n + 1
		neighbor = xs[j]
		s = interior[j]

		foo = first(r) - 2overlap - stencil:first(r) - 2overlap - 1
		bar = last(s) - stencil + 1:last(s)

		for (k, l) in zip(foo, bar)
			center[k] = neighbor[l]
		end

		# right
		j = (i) % n + 1
		neighbor = xs[j]
		s = interior[j]

		foo = last(r) + 2overlap + 1:last(r) + 2overlap + stencil
		bar = first(s):first(s) + stencil - 1

		for (k, l) in zip(foo, bar)
			center[k] = neighbor[l]
		end
	end

	nothing
end

# ╔═╡ a8f7eeb7-799a-48bc-841e-9d35e99b0268
update!(bs, b)

# ╔═╡ 8302e1ed-d953-43b9-afbb-4e45a5cf628b
update!(Interior, bs, b)

# ╔═╡ 5803f096-760c-4ec3-934e-72f74bc72308
bs[begin]

# ╔═╡ fa532477-fec7-47de-909d-f07713eccb7b
parent(bs[begin + 1])

# ╔═╡ 51688333-7f36-49a3-809e-0cf7338b6f9c
begin
	n = size(bs, 1)

	i = 2
	local center = bs[i]
	r = interior[i]

	local overlap = getoverlap(center)
	local stencil = getstencil(center)

	j = (i) % n + 1
	neighbor = bs[j]
	s = interior[j]
	@show j
	#=

	foo = first(r) - 2overlap - stencil:first(r) - 2overlap - 1
	local bar = last(s) - stencil + 1:last(s)

	@show foo, bar
	=#
#	center[foo] .= neighbor[bar]
#=	for (k, l) in zip(foo, bar)
		center[k] = neighbor[l]
	end
=#
end

# ╔═╡ ba8f8824-647e-4897-b551-6b6a0e3244e3
interior[1], interior[2]

# ╔═╡ f9590a9f-74d1-49a7-9c39-db2a9f28a09c


# ╔═╡ 9bb8013c-c84a-448c-bc8f-7727c7d6a6db
axes(bs[1])

# ╔═╡ 8438a49c-73a0-4b76-b33b-d8b180a4c478
axes(parent(bs[2]))

# ╔═╡ 3882cba8-d9fa-4b07-83cc-1596a4796bcc
cycle

# ╔═╡ 602bd53f-186c-41a5-b9a2-7d351d8544fb
begin
	local i = 9
	local n = 9
	(n + i - 2) % n + 1, (i) % n + 1
end

# ╔═╡ f811ca17-7631-4c1b-80c6-2fb3e15b5d83


# ╔═╡ 0c908d33-5253-4511-8536-abf1c9d696fa
update!(Ghost, bs, interior)

# ╔═╡ 6753649a-a772-4f2b-bd90-1e5d866664b9
begin
	local n = 16
	x = [sinpi(i / n) for i in 1:n]
	y = x
end

# ╔═╡ 917a4a9e-e8d3-4815-88c0-4e2435b89435
x

# ╔═╡ 561004ca-a30c-4f26-a024-ca54afc3f997
y[2] = -1

# ╔═╡ 351fb0f1-8b36-473c-a305-9aa97caaaafa
xs = [x; x; x]

# ╔═╡ b080fc54-9526-4ff2-8a18-5f458c288bce
vcat(x, xs, x)

# ╔═╡ 63ed6044-5365-40b4-8f3d-1cb442dae4d1
temp = CartesianIndices((2:64,))

# ╔═╡ Cell order:
# ╠═110d3d54-9225-4a44-83ce-3a7f509232fd
# ╟─06626098-d059-11eb-0434-6f829bce295f
# ╠═e6a47c50-0d83-4333-af8e-6b493aef773c
# ╠═74880f2a-cd4b-49af-a574-c680785b9be6
# ╠═8683a8d9-af47-4c3b-a075-b7b4796c297f
# ╠═5429496f-9473-48d2-9842-50bbace1403c
# ╠═e630dee3-1781-4d81-9185-45096a87047c
# ╠═2b2cf2c0-c38d-4991-8a33-7a3b481754ec
# ╠═087bc9b7-817d-49fc-ab00-513ab03f342a
# ╠═f30e6e06-ca06-4a0e-8dc8-a9b45eb6db5c
# ╟─a8d9eeff-bd77-4e02-89b0-3f676fa11fc3
# ╠═49bf2458-e192-4c01-80e4-f48b37bca3d2
# ╠═cb9a7522-8dba-4c95-a7c3-93ee8f2fd21b
# ╠═baaa29dd-1fb6-451b-961c-20c9a45be294
# ╠═d8452319-d87b-4d11-b21f-f2c8306179ee
# ╠═c72529da-9395-4409-9f55-b15e46611e0d
# ╠═a1a7383b-9585-4952-9a27-a6ed74b60ce0
# ╠═3ecd7b07-e603-4994-8941-b704bae05299
# ╠═0d2d52a4-c656-41dd-b641-17c209c3affa
# ╠═2d8bf534-0203-4d58-9f2d-1aabd99176a3
# ╠═0fecb93e-8d74-4c9d-a5b6-7fe43e175268
# ╠═449cd7fb-77f7-4c53-b665-0d322aed465d
# ╠═d8298fed-3266-47c7-8b0a-1361a940c557
# ╠═a615ce0f-f57e-4e50-9334-4b487f9a2443
# ╠═37eb5392-9c17-4d06-883b-0d3693a954fc
# ╠═b39d2736-c827-452b-a5a3-adbf19eb781d
# ╠═3c59a707-3b63-450d-90e1-1ed6c00da6c2
# ╠═86247385-1aed-4ed6-b70f-4a567201889e
# ╠═cd8534b0-e522-4d42-8cc0-bf3ae75ad92e
# ╠═47ce9c5c-4a82-43dc-b529-d83e55532e4c
# ╠═0eba4d39-e3f7-45a1-96d6-847109e48e23
# ╠═1f31f525-afc4-4e20-ae3e-fd3fb93aca24
# ╠═46277571-c698-4dd7-a6a7-a6c609d0845e
# ╠═5dab8c53-4fac-4ac4-833c-78918dd3d823
# ╠═5da69b07-3d42-4ee1-86ad-debee168de8d
# ╠═1542f9dd-21e7-4e4e-a50e-5a234b1f8abe
# ╠═dc293c44-578f-4d45-aca0-a6f766818314
# ╠═f5d9153c-e31d-4562-a05b-19b4b7f0ca01
# ╠═a8f7eeb7-799a-48bc-841e-9d35e99b0268
# ╠═5182c5ca-b1f0-4773-97fe-33e9a939c540
# ╠═49118046-3d81-4e2c-80cf-0e438e3e391a
# ╠═a00e640f-cec0-4a17-8215-a5d59e42fe07
# ╠═798003a3-a0bc-47b8-b4ad-77a112da466e
# ╠═8302e1ed-d953-43b9-afbb-4e45a5cf628b
# ╠═675b43e1-66ac-4143-8b9a-ae2ff3bad1d0
# ╠═5803f096-760c-4ec3-934e-72f74bc72308
# ╠═fa532477-fec7-47de-909d-f07713eccb7b
# ╠═51688333-7f36-49a3-809e-0cf7338b6f9c
# ╠═ba8f8824-647e-4897-b551-6b6a0e3244e3
# ╠═f9590a9f-74d1-49a7-9c39-db2a9f28a09c
# ╠═9bb8013c-c84a-448c-bc8f-7727c7d6a6db
# ╠═8438a49c-73a0-4b76-b33b-d8b180a4c478
# ╠═3882cba8-d9fa-4b07-83cc-1596a4796bcc
# ╠═602bd53f-186c-41a5-b9a2-7d351d8544fb
# ╠═f811ca17-7631-4c1b-80c6-2fb3e15b5d83
# ╠═0c908d33-5253-4511-8536-abf1c9d696fa
# ╠═6753649a-a772-4f2b-bd90-1e5d866664b9
# ╠═917a4a9e-e8d3-4815-88c0-4e2435b89435
# ╠═561004ca-a30c-4f26-a024-ca54afc3f997
# ╠═351fb0f1-8b36-473c-a305-9aa97caaaafa
# ╠═b080fc54-9526-4ff2-8a18-5f458c288bce
# ╠═63ed6044-5365-40b4-8f3d-1cb442dae4d1
