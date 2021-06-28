### A Pluto.jl notebook ###
# v0.14.2

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

!!! note "Note"

	Could this be generalized to non-unit spacing? Also, `IdOffsetRange` can be a bit of a pain...


"""

# ╔═╡ e6a47c50-0d83-4333-af8e-6b493aef773c
struct UniformPartition{N,T,R<:NTuple{N,AbstractUnitRange{T}}} <: AbstractArray{R,N}
	indices::R
	nover::NTuple{N,T}
	nproc::NTuple{N,T}
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

	ls = first.(indices)
	hs = last.(indices)
	Ts = typeof.(indices)

	map(Ts, ls, hs, nover, nproc, index) do T, l, h, o, p, i
		d = h - l + 1
		T(
			max(d * (i - 1) ÷ p + l - o ÷ 2, l),
			min(d * i ÷ p - 1 + l + (o + 1) ÷ 2, h)
		)
	end
end

# ╔═╡ a8d9eeff-bd77-4e02-89b0-3f676fa11fc3
md"""
# Padding

"""

# ╔═╡ a9d0b495-e093-45d0-8f11-60ec2b8f6f0d
function pad(indices, n)
	UnitRange(
		first(indices) - n,
		last(indices) + n
	)
end

# ╔═╡ cb9a7522-8dba-4c95-a7c3-93ee8f2fd21b
md"""
# Ghosting

"""

# ╔═╡ baaa29dd-1fb6-451b-961c-20c9a45be294
struct GhostedArray{T,N,AA<:AbstractArray{T,N}} <: AbstractArray{T,N}
	parent::AA
	nghost::NTuple{N,Int}
end

# ╔═╡ a70d400a-fe4a-46a8-a4df-c5b9abe78c4d
function ghostedarray(::Type{T}, indices, n) where {T}
	ghosted = pad.(indices, n)
	data = OffsetArray(Array{T}(undef, length.(ghosted)), ghosted)
	GhostedArray(data, n)
end

# ╔═╡ d8452319-d87b-4d11-b21f-f2c8306179ee
Base.parent(a::GhostedArray) = a.parent

# ╔═╡ c09cc1ad-eeba-4be0-850e-b28d846b478b
getghost(a::GhostedArray) = a.nghost

# ╔═╡ c72529da-9395-4409-9f55-b15e46611e0d
Base.getindex(a::GhostedArray, i...) = getindex(parent(a), i...)

# ╔═╡ a1a7383b-9585-4952-9a27-a6ed74b60ce0
Base.setindex!(a::GhostedArray, i...) = setindex!(parent(a), i...)

# ╔═╡ 3ecd7b07-e603-4994-8941-b704bae05299
function Base.axes(a::GhostedArray)
	nghost = getghost(a)
	indices = axes(parent(a))
	pad.(indices, .-nghost)
end

# ╔═╡ 0d2d52a4-c656-41dd-b641-17c209c3affa
Base.size(a::GhostedArray) = length.(axes(a))

# ╔═╡ 00a88eb8-2d1e-4915-b6bf-22fb12fc801e
md"""
# Communication

"""

# ╔═╡ b564ec93-b776-4a71-a2e5-0361c6178d8d
function update!(ys, x)
	for y in ys
		z = parent(y)
		indices = CartesianIndices(z)
		z .= x[indices]
	end
end

# ╔═╡ 8ecb2268-e0db-4edf-bae2-5cbe552b5e25
md"""
# Example

Let's first create a partition.

"""

# ╔═╡ f30e6e06-ca06-4a0e-8dc8-a9b45eb6db5c
partition = UniformPartition(
	# indices
	(32:64, 17:32),
	# overlap (nover)
	(1, 2),
	# processes (nproc)
	(9, 3)
)

# ╔═╡ 37bfdb62-a381-49f2-8ed2-0b8997bcf25b
md"""
For second order finite difference, a single ghost node suffices.

"""

# ╔═╡ 2d8bf534-0203-4d58-9f2d-1aabd99176a3
nghost = (1, 1)

# ╔═╡ 2d6fad34-e6c2-490c-9275-d0a0e74c8717
md"""
We can now create global and local arrays.

"""

# ╔═╡ 0fecb93e-8d74-4c9d-a5b6-7fe43e175268
b = ghostedarray(Float64, getindices(partition), nghost)

# ╔═╡ e6d64754-2e1a-4567-9b0f-5229b4b3eb4b
bs = map(partition) do indices
	ghostedarray(Float64, indices, nghost)
end

# ╔═╡ 3ab3345d-f318-4052-a3ac-4254e5f37dca
md"""
And transfer information from global to local...

"""

# ╔═╡ ffdb2781-5a40-495c-a73b-737b98f5a79d
update!(bs, b)

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
# ╟─a8d9eeff-bd77-4e02-89b0-3f676fa11fc3
# ╠═a9d0b495-e093-45d0-8f11-60ec2b8f6f0d
# ╟─cb9a7522-8dba-4c95-a7c3-93ee8f2fd21b
# ╠═baaa29dd-1fb6-451b-961c-20c9a45be294
# ╠═a70d400a-fe4a-46a8-a4df-c5b9abe78c4d
# ╠═d8452319-d87b-4d11-b21f-f2c8306179ee
# ╠═c09cc1ad-eeba-4be0-850e-b28d846b478b
# ╠═c72529da-9395-4409-9f55-b15e46611e0d
# ╠═a1a7383b-9585-4952-9a27-a6ed74b60ce0
# ╠═3ecd7b07-e603-4994-8941-b704bae05299
# ╠═0d2d52a4-c656-41dd-b641-17c209c3affa
# ╟─00a88eb8-2d1e-4915-b6bf-22fb12fc801e
# ╠═b564ec93-b776-4a71-a2e5-0361c6178d8d
# ╟─8ecb2268-e0db-4edf-bae2-5cbe552b5e25
# ╠═f30e6e06-ca06-4a0e-8dc8-a9b45eb6db5c
# ╟─37bfdb62-a381-49f2-8ed2-0b8997bcf25b
# ╠═2d8bf534-0203-4d58-9f2d-1aabd99176a3
# ╟─2d6fad34-e6c2-490c-9275-d0a0e74c8717
# ╠═0fecb93e-8d74-4c9d-a5b6-7fe43e175268
# ╠═e6d64754-2e1a-4567-9b0f-5229b4b3eb4b
# ╟─3ab3345d-f318-4052-a3ac-4254e5f37dca
# ╠═ffdb2781-5a40-495c-a73b-737b98f5a79d
