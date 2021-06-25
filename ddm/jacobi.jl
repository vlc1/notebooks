### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 110d3d54-9225-4a44-83ce-3a7f509232fd
using OffsetArrays

# ╔═╡ 1d7ff3fd-e44e-4f13-ba6b-e6902dbb2fbb
using Plots

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
	struct UniformPartition{N,T} <: AbstractArray{CartesianIndices{N,T},N}
		ndof::NTuple{N,Int} # also represent this with CartesianIndices
		nproc::NTuple{N,Int}
	end

	UniformPartition(ndof::T, nproc::T) where {N,T<:NTuple{N}} =
		UniformPartition{N,NTuple{N,UnitRange{Int}}}(ndof, nproc)
end

# ╔═╡ 74880f2a-cd4b-49af-a574-c680785b9be6
getdof(p::UniformPartition) = p.ndof

# ╔═╡ 5429496f-9473-48d2-9842-50bbace1403c
getproc(p::UniformPartition) = p.nproc

# ╔═╡ e630dee3-1781-4d81-9185-45096a87047c
Base.size(p::UniformPartition) = getproc(p)

# ╔═╡ 2b2cf2c0-c38d-4991-8a33-7a3b481754ec
Base.getindex(part::UniformPartition, index::CartesianIndex) =
	getindex(part, index.I...)

# ╔═╡ 087bc9b7-817d-49fc-ab00-513ab03f342a
function Base.getindex(part::UniformPartition, index...)
	ndof, nproc = getdof(part), getproc(part)

	map(ndof, nproc, index) do d, p, i
		UnitRange(
			floor(Int, d * (i - 1) / p) + 1,
			floor(Int, d * i / p) + 1 - 1
		)
	end |> CartesianIndices
end

# ╔═╡ e92f485c-e2a4-45b8-8f56-8afa3aff4b08
partition = UniformPartition((64,32), (9,5))

# ╔═╡ cb9a7522-8dba-4c95-a7c3-93ee8f2fd21b
md"""
# Domains

"""

# ╔═╡ 2d8bf534-0203-4d58-9f2d-1aabd99176a3
begin
	local nover, nst = 1, 1
	param = (nover, nst)
end

# ╔═╡ 429fe85b-87ac-4ad5-9512-df7130d3bf48
abstract type AbstractDomain end

# ╔═╡ c2b83803-de47-405c-bc35-a8318c4d538d
md"""
## Interior domain

"""

# ╔═╡ 97603474-d72a-4841-a206-09bd1242d608
struct Interior <: AbstractDomain end

# ╔═╡ 105935ef-8d8e-4934-8485-9e7a41549c34
getshift(::Interior, nover, args...) = -nover

# ╔═╡ 9e192127-c4ef-4924-a2af-976fe7ac52a0
interior(::AbstractArray{T,N}) where {T,N} =
	ntuple(i -> Interior(), Val(N))

# ╔═╡ 7a05de82-ef35-493f-9183-350c5582b58d
md"""
## Half overlap

"""

# ╔═╡ 81860fb4-ba52-4ccc-90cd-318f2124683f
struct SemiOverlap <: AbstractDomain end

# ╔═╡ 68656ff6-4632-4193-b1f9-a026c8ae1a0e
getshift(::SemiOverlap, args...) = 0

# ╔═╡ 6dfa3945-59d5-4425-8469-e5d513cf7133
semioverlap(::AbstractArray{T,N}) where {T,N} =
	ntuple(i -> SemiOverlap(), Val(N))

# ╔═╡ c8327bb8-5f3f-46b5-9a5e-d489f53eeac8
md"""
## Full overlap

"""

# ╔═╡ 859faf02-ef8f-4a25-95d5-3d7695ec64f9
struct FullOverlap <: AbstractDomain end

# ╔═╡ d0a4c4e4-ac73-4e5a-a6b1-f188c0e5de1e
fulloverlap(::AbstractArray{T,N}) where {T,N} =
	ntuple(i -> FullOverlap(), Val(N))

# ╔═╡ 893975d5-c3e3-4289-94a8-62f9fb7d0c3c
getshift(::FullOverlap, nover, args...) = nover

# ╔═╡ adce8693-a881-47d3-9f31-b68ff9c3b617
md"""
## Ghosted overlap

"""

# ╔═╡ 6b5cb70a-6b1d-42d8-8eb1-d6d63da4a347
struct GhostedOverlap <: AbstractDomain end

# ╔═╡ e91c9a52-1a27-40b2-b822-14f99d293894
getshift(::GhostedOverlap, nover, nst) = nover + nst

# ╔═╡ 9f56bc48-43bc-4fce-a9b3-a2eed86b6ccd
ghostedoverlap(::AbstractArray{T,N}) where {T,N} =
	ntuple(i -> GhostedOverlap(), Val(N))

# ╔═╡ c172fbca-737d-4c4c-a75a-3ef9444f3d48
md"""
## Selector

"""

# ╔═╡ 2a3e3384-2479-4e9e-a72d-7e2b88d5de49
function select(domains, indices, args...)
	map(domains, Tuple(first(indices)), Tuple(last(indices))) do domain, start, stop
		shift = getshift(domain, args...)
		start - shift:stop + shift
	end |> CartesianIndices
end

# ╔═╡ 09b19bb8-05c0-4373-bb10-81b6388cc872
function interior(indices, args...)
	domains = interior(indices)
	select(domains, indices, args...)
end

# ╔═╡ 852d21fa-55a6-4541-bea0-f7cf7c0eee37
function semioverlap(indices, args...)
	domains = semioverlap(indices)
	select(domains, indices, args...)
end

# ╔═╡ b4b0114a-0cb5-458b-bc08-646b07c447b4
function fulloverlap(indices, args...)
	domains = fulloverlap(indices)
	select(domains, indices, args...)
end

# ╔═╡ 2e0c9b29-609f-467a-b0d6-8f57afcb3b5a
function ghostedoverlap(indices, args...)
	domains = ghostedoverlap(indices)
	select(domains, indices, args...)
end

# ╔═╡ 7299bb9c-cf6e-4840-a3f6-e3ed440d78ed
md"""
To get the interior of a partition :

"""

# ╔═╡ a776f655-5c0b-4e55-be59-7575f10bd5d2
interior(partition[2,3], param...)

# ╔═╡ ca7660d4-1faa-4fe6-927f-408cdd415c5f
md"""
# Global and local arrays

*À la* PETSc...

"""

# ╔═╡ aee7038c-3d40-45aa-ae54-9333ecffeeb4
begin
	local n = getdof(partition)
	local x = [[(2i - 1) / n[j] for i in 1:n[j]] for j in 1:length(n)]
	b = reshape(
		[sinpi(x[1][i]) * cospi(x[2][j]) for i in 1:n[1] for j in 1:n[2]],
		n
		)
	#scatter(x[1], b, label = "rhs")
	#x
end

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

# ╔═╡ 7510547d-bfc9-4b9e-918d-892251f0c14a
#Base.parent(a::DDMArray) = a.parent

# ╔═╡ 0eba4d39-e3f7-45a1-96d6-847109e48e23
#overlap(a::DDMArray) = a.nover

# ╔═╡ 1f31f525-afc4-4e20-ae3e-fd3fb93aca24
#stencil(a::DDMArray) = a.nst

# ╔═╡ 93403c71-e61e-4a2e-9277-ff51fb009637
#Base.getindex(a::DDMArray, i...) = getindex(parent(a), i...)

# ╔═╡ 00e0e5e0-bfe1-41d6-8950-f5a52af1ad66
#Base.setindex!(a::DDMArray, i...) = setindex!(parent(a), i...)

# ╔═╡ ae492709-aba1-473f-8467-80a5146a7be3
#Base.axes(a::DDMArray) = axes(parent(a))

# ╔═╡ c430985d-3bef-47ff-a665-1850e0239231
#Base.size(a::DDMArray) = size(parent(a))

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
function initialize(x, partition, param...)
	map(partition) do indices
		ghosted = ghostedoverlap(indices, param...)
		OffsetArray(similar(x, axes(ghosted)), ghosted)
	end
end

# ╔═╡ dc293c44-578f-4d45-aca0-a6f766818314
bs = initialize(b, partition, param...)

# ╔═╡ f5d9153c-e31d-4562-a05b-19b4b7f0ca01
function update!(bs, b, param...)
	
end

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

# ╔═╡ Cell order:
# ╠═110d3d54-9225-4a44-83ce-3a7f509232fd
# ╟─06626098-d059-11eb-0434-6f829bce295f
# ╠═e6a47c50-0d83-4333-af8e-6b493aef773c
# ╠═74880f2a-cd4b-49af-a574-c680785b9be6
# ╠═5429496f-9473-48d2-9842-50bbace1403c
# ╠═e630dee3-1781-4d81-9185-45096a87047c
# ╠═2b2cf2c0-c38d-4991-8a33-7a3b481754ec
# ╠═087bc9b7-817d-49fc-ab00-513ab03f342a
# ╠═e92f485c-e2a4-45b8-8f56-8afa3aff4b08
# ╟─cb9a7522-8dba-4c95-a7c3-93ee8f2fd21b
# ╠═2d8bf534-0203-4d58-9f2d-1aabd99176a3
# ╠═429fe85b-87ac-4ad5-9512-df7130d3bf48
# ╟─c2b83803-de47-405c-bc35-a8318c4d538d
# ╠═97603474-d72a-4841-a206-09bd1242d608
# ╠═105935ef-8d8e-4934-8485-9e7a41549c34
# ╠═9e192127-c4ef-4924-a2af-976fe7ac52a0
# ╠═09b19bb8-05c0-4373-bb10-81b6388cc872
# ╟─7a05de82-ef35-493f-9183-350c5582b58d
# ╠═81860fb4-ba52-4ccc-90cd-318f2124683f
# ╠═68656ff6-4632-4193-b1f9-a026c8ae1a0e
# ╠═6dfa3945-59d5-4425-8469-e5d513cf7133
# ╠═852d21fa-55a6-4541-bea0-f7cf7c0eee37
# ╟─c8327bb8-5f3f-46b5-9a5e-d489f53eeac8
# ╠═859faf02-ef8f-4a25-95d5-3d7695ec64f9
# ╠═d0a4c4e4-ac73-4e5a-a6b1-f188c0e5de1e
# ╠═893975d5-c3e3-4289-94a8-62f9fb7d0c3c
# ╠═b4b0114a-0cb5-458b-bc08-646b07c447b4
# ╟─adce8693-a881-47d3-9f31-b68ff9c3b617
# ╠═6b5cb70a-6b1d-42d8-8eb1-d6d63da4a347
# ╠═e91c9a52-1a27-40b2-b822-14f99d293894
# ╠═9f56bc48-43bc-4fce-a9b3-a2eed86b6ccd
# ╠═2e0c9b29-609f-467a-b0d6-8f57afcb3b5a
# ╟─c172fbca-737d-4c4c-a75a-3ef9444f3d48
# ╠═2a3e3384-2479-4e9e-a72d-7e2b88d5de49
# ╟─7299bb9c-cf6e-4840-a3f6-e3ed440d78ed
# ╠═a776f655-5c0b-4e55-be59-7575f10bd5d2
# ╟─ca7660d4-1faa-4fe6-927f-408cdd415c5f
# ╠═1d7ff3fd-e44e-4f13-ba6b-e6902dbb2fbb
# ╠═aee7038c-3d40-45aa-ae54-9333ecffeeb4
# ╠═3c59a707-3b63-450d-90e1-1ed6c00da6c2
# ╠═86247385-1aed-4ed6-b70f-4a567201889e
# ╠═cd8534b0-e522-4d42-8cc0-bf3ae75ad92e
# ╠═47ce9c5c-4a82-43dc-b529-d83e55532e4c
# ╠═7510547d-bfc9-4b9e-918d-892251f0c14a
# ╠═0eba4d39-e3f7-45a1-96d6-847109e48e23
# ╠═1f31f525-afc4-4e20-ae3e-fd3fb93aca24
# ╠═93403c71-e61e-4a2e-9277-ff51fb009637
# ╠═00e0e5e0-bfe1-41d6-8950-f5a52af1ad66
# ╠═ae492709-aba1-473f-8467-80a5146a7be3
# ╠═c430985d-3bef-47ff-a665-1850e0239231
# ╠═46277571-c698-4dd7-a6a7-a6c609d0845e
# ╠═5dab8c53-4fac-4ac4-833c-78918dd3d823
# ╠═dc293c44-578f-4d45-aca0-a6f766818314
# ╠═f5d9153c-e31d-4562-a05b-19b4b7f0ca01
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
