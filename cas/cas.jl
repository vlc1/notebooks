### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 09455e66-a91d-11eb-15fa-d914f4279d8f
using SymbolicUtils

# ╔═╡ 20a40aa9-0a3f-4bbf-9338-9a8fa6704dcd
using Symbolics

# ╔═╡ 56d01643-bf98-4e92-b439-a0ce2b693fb7
using LinearAlgebra

# ╔═╡ a53bcfac-78c4-4270-8a74-45edf76bf1db
md"""
# Definitions

"""

# ╔═╡ d6d4c405-3266-4dde-8930-f16ba9aa7811
#@syms λ::Number μ::Number

# ╔═╡ cbf1fcfc-8a85-49fe-aef5-cd16dff3e735
@syms u::AbstractVector v::AbstractVector w::AbstractVector

# ╔═╡ 13c973de-6b3c-44e5-97bc-a8f0ecf47ed0
#@syms (*)(::Number, ::AbstractVector)::AbstractVector

# ╔═╡ 9f344936-cf57-4d98-b769-75ff3a848d4b
#λ * u

# ╔═╡ cd055cb0-a0cb-4293-91ba-b004ef59363c
# Not the right way to do it
@syms (+)(::Vararg{AbstractVector})::AbstractVector

# ╔═╡ 2af7c994-6cb3-4aca-add5-c3f207cc42ee
@syms (*)(::Vararg{AbstractVector})::AbstractVector

# ╔═╡ db197bf0-b9c4-499f-bae3-adb268636527
u + v + w

# ╔═╡ 004e6dc6-f175-49bd-b53c-dd3daf11de36
#λ * u

# ╔═╡ d99b9192-b846-4034-a708-9a37213a4de0
@syms δ(x::AbstractVector)::AbstractVector

# ╔═╡ 55e249c2-5f5d-4f7d-94ab-7d8b9105566c
@syms σ(x::AbstractVector)::AbstractVector

# ╔═╡ 18ee34e3-293f-4e76-8598-f03048b33061
σ(δ(u))

# ╔═╡ eeb514f0-5b1a-4b88-9be2-2df77325c1d4
isequal(u + v, v + u)

# ╔═╡ e656aaf3-3679-495d-b734-da1e0df40519
md"""
# Rules

## Linearity

"""

# ╔═╡ 81ec6e5b-36d0-4696-b3c5-45235c80818d
r1 = @rule(σ(+(~~us)) => reduce(+, map(σ, ~~us)))

# ╔═╡ cba26072-2db1-4ca4-9565-f21f7b7befb8
r1(σ(u + v))

# ╔═╡ 7821162a-3dde-4f05-9d39-f49d4920d641
r1(σ(u + v + w))

# ╔═╡ 11440358-dd84-4c9f-bd28-9a9878620a40
r2 = @rule(δ(+(~~us)) => reduce(+, map(δ, ~~us)))

# ╔═╡ 5e26e55e-b876-4949-afe0-3d41e6414cab
r2(δ(u + v))

# ╔═╡ e90dc1e3-63a6-456e-b4dc-8b7eaf37c20b
r2(δ(u + v + w))

# ╔═╡ 4d323dbc-87a5-4921-9861-0401f500dda8
#r3 = @rule(σ((~λ) * (~u)) => (~λ) * σ(~u))

# ╔═╡ 7685fdb2-b7da-48b3-ae96-94e7ceaafb5b
#r3(σ(λ * u))

# ╔═╡ 4583d81f-58b2-44fa-b3be-ece2a09ae475
#r4 = @rule(δ((~λ) * (~u)) => (~λ) * δ(~u))

# ╔═╡ 75f1814e-97fe-4a21-a9ec-cf9198ae7993
#r4(δ(λ * u))

# ╔═╡ 26efc1e3-5ab7-46b6-91ad-37da97c6d8d0
@rule 

# ╔═╡ 45a9227c-ac8f-4b54-a55d-7c448984a1cc
u * v

# ╔═╡ 6e2598b2-83e8-4df8-9ef1-08a869b4cb6a
#@syms (.*)(::Vararg{AbstractVector})::AbstractVector

# ╔═╡ 34dc62ea-0038-4da6-a922-029bc3e6492b
md"""
# Commutativity

"""

# ╔═╡ a2de2dc3-5ad2-4e8e-8182-cb1e9111dd5f
r1(σ(u + v))

# ╔═╡ 65c9ea39-03ae-475d-9ecc-94fd3f372e5e
r2(δ(u + v))

# ╔═╡ 1170f0c9-f8a8-4319-a75a-38559cf7f91c
#r3(δ(λ * (u + v)))

# ╔═╡ 81aebcd9-38ab-445e-b4f1-72b65a52f91b
md"""
# Kronecker product

"""

# ╔═╡ 366a8f9a-3452-4b77-a406-50995378a899
@syms I::UniformScaling

# ╔═╡ 53122d4a-0610-458b-b04f-357476d20e32
@syms kron(x, y)

# ╔═╡ 6ae06407-e55d-4c65-a4b0-6f99eed75eb1
md"""
# Adjoint

"""

# ╔═╡ 07a1c71d-3038-44eb-949d-ee0728c21014
@syms δ̲(::AbstractVector)::AbstractVector

# ╔═╡ Cell order:
# ╠═09455e66-a91d-11eb-15fa-d914f4279d8f
# ╠═20a40aa9-0a3f-4bbf-9338-9a8fa6704dcd
# ╟─a53bcfac-78c4-4270-8a74-45edf76bf1db
# ╠═d6d4c405-3266-4dde-8930-f16ba9aa7811
# ╠═cbf1fcfc-8a85-49fe-aef5-cd16dff3e735
# ╠═13c973de-6b3c-44e5-97bc-a8f0ecf47ed0
# ╠═9f344936-cf57-4d98-b769-75ff3a848d4b
# ╠═cd055cb0-a0cb-4293-91ba-b004ef59363c
# ╠═2af7c994-6cb3-4aca-add5-c3f207cc42ee
# ╠═db197bf0-b9c4-499f-bae3-adb268636527
# ╠═004e6dc6-f175-49bd-b53c-dd3daf11de36
# ╠═d99b9192-b846-4034-a708-9a37213a4de0
# ╠═55e249c2-5f5d-4f7d-94ab-7d8b9105566c
# ╠═18ee34e3-293f-4e76-8598-f03048b33061
# ╠═eeb514f0-5b1a-4b88-9be2-2df77325c1d4
# ╟─e656aaf3-3679-495d-b734-da1e0df40519
# ╠═81ec6e5b-36d0-4696-b3c5-45235c80818d
# ╠═cba26072-2db1-4ca4-9565-f21f7b7befb8
# ╠═7821162a-3dde-4f05-9d39-f49d4920d641
# ╠═11440358-dd84-4c9f-bd28-9a9878620a40
# ╠═5e26e55e-b876-4949-afe0-3d41e6414cab
# ╠═e90dc1e3-63a6-456e-b4dc-8b7eaf37c20b
# ╠═4d323dbc-87a5-4921-9861-0401f500dda8
# ╠═7685fdb2-b7da-48b3-ae96-94e7ceaafb5b
# ╠═4583d81f-58b2-44fa-b3be-ece2a09ae475
# ╠═75f1814e-97fe-4a21-a9ec-cf9198ae7993
# ╠═26efc1e3-5ab7-46b6-91ad-37da97c6d8d0
# ╠═45a9227c-ac8f-4b54-a55d-7c448984a1cc
# ╠═6e2598b2-83e8-4df8-9ef1-08a869b4cb6a
# ╟─34dc62ea-0038-4da6-a922-029bc3e6492b
# ╠═a2de2dc3-5ad2-4e8e-8182-cb1e9111dd5f
# ╠═65c9ea39-03ae-475d-9ecc-94fd3f372e5e
# ╠═1170f0c9-f8a8-4319-a75a-38559cf7f91c
# ╠═81aebcd9-38ab-445e-b4f1-72b65a52f91b
# ╠═56d01643-bf98-4e92-b439-a0ce2b693fb7
# ╠═366a8f9a-3452-4b77-a406-50995378a899
# ╠═53122d4a-0610-458b-b04f-357476d20e32
# ╠═6ae06407-e55d-4c65-a4b0-6f99eed75eb1
# ╠═07a1c71d-3038-44eb-949d-ee0728c21014
