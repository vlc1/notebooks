### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ c810af14-b553-4a25-9423-c18850dae690
using IterativeSolvers

# ╔═╡ 626f475f-bc14-4f82-bac8-3976bf3f25d0
using DistributedArrays

# ╔═╡ 42bbbcd8-1d50-42df-a700-42743d10e417
using LinearAlgebra

# ╔═╡ 6dec645a-8a0b-49c8-8ae9-37c7c7f2b25e
md"""
# Matrices diagonales

Soit une matrice diagonale, ``M \left ( x \right )`` de diagonale ``x \in \mathbb{R} ^ n``.

Le module `LinearAlgebra` de la bibliothèque standard a son propre type (`Diagonal`) pour représenter ces matrices, mais on définira notre propre type ci-dessous, `MyScaling`.

"""

# ╔═╡ c0340c58-c233-11eb-32fd-ddb8506d9ecb
begin
	struct MyDiagonal{T, S<:AbstractVector} <: AbstractMatrix{T}
		diag::S
	end
	MyDiagonal(x) = MyDiagonal{eltype(x), typeof(x)}(x)
end

# ╔═╡ 157ab7ce-ed31-430e-a230-438f02630753
md"""
L'affichage d'objets de type `MyScaling` demande *a minima* la définition des méthodes suivantes.

"""

# ╔═╡ da7cad0e-40eb-404f-afdf-f3e1f5aa2307
Base.size(A::MyDiagonal) = ntuple(i -> size(A.diag, 1), 2)

# ╔═╡ 23e036d6-c67d-4e11-9d17-b87f6f0657c8
Base.getindex(A::MyDiagonal, i, j) = i == j ? A.diag[i] : zero(eltype(A))

# ╔═╡ 5029cc06-922e-4c74-953d-cc660bb551c2
md"""
La matrice
```math
\left ( \begin{matrix}
\sqrt{2} & 0 & 0 \\
0 & -1 & 0 \\
0 & 0 & 1 / 3
\end{matrix} \right ) 
```
peut alors être construite comme suit.

"""

# ╔═╡ 9fd097d5-1696-4426-a2a0-e78e7d611df5
foo = MyDiagonal([√2, -1, 1/3])

# ╔═╡ a5cf0863-b6a7-4dab-ae80-68f58484bf17
md"""
La multiplication se fait tout simplement en ajoutant la méthode suivante à la fonction `*` du module `Base`.

"""

# ╔═╡ 90f90d71-452f-45e0-bab9-3443b8e9cd98
Base.:*(A::MyDiagonal, x::AbstractVector) = A.diag .* x

# ╔═╡ efac4950-d829-45b2-a9dc-958fb9d57996
md"""
Ça marche comme prévu.

"""

# ╔═╡ d46458b9-0b34-4fb7-bf94-442db95b9de1
bar = rand(3)

# ╔═╡ 159d4b51-50e3-4955-833a-1510b321ec61
foo * bar

# ╔═╡ dbb335aa-b9b2-40e8-a431-a05a5639ea2b
n = 8

# ╔═╡ e5ee30dd-633a-4c37-bf78-92a1a21d65cd
md"""
# Le *package* `IterativeSolvers.jl`

Voyons voir si ça marche avec le CG de `IterativeSolvers.jl`. Soit ``A`` une matrice diagonal positive definie de taille $(n).

"""

# ╔═╡ a61234f3-72d1-4d7a-b870-ef89542e7567
A = MyDiagonal(rand(n) .^ 2 .+ √eps())

# ╔═╡ dc49be41-b4b5-4fee-8e7f-c45969728e71
b = rand(n)

# ╔═╡ eeaf4a3e-51f5-4951-a0a5-0f1e6c16c2cd
x = cg(A, b)

# ╔═╡ 38c3e98c-4cbc-4559-982c-6b9f94574ad2
if isapprox(b, A * x; atol = ∛eps())
	md"""
	!!! tip "😃 ça marche !"
	"""
else
	md"""
	!!! danger "😠 ça marche pas !"
	"""
end

# ╔═╡ aa75ad03-74d9-4602-9219-e641f28a3a15
m = 100

# ╔═╡ b1b12b1b-1cc4-4c1d-9352-95551d1b0aa3
md"""
# Le *package* `DistributedArrays.jl`

Voyons voir maintenant si on peut refaire la même chose, mais en distribuant des vecteurs de taille $m.

"""

# ╔═╡ 64b6cc94-3297-4a89-8c6f-23f9da1e110d
dA = MyDiagonal(drand(m))

# ╔═╡ cb3eb1f1-c1e3-4b08-8243-e8d523cc355c
db = drand(m)

# ╔═╡ b8f68ab8-84eb-4a11-bb5c-35d8b5f95057
dx = cg(dA, db)

# ╔═╡ f61420f4-024e-4bcf-99dc-881ccccc5404
md"""
# Avec la bibliothèque standard

Ça ne marche pas non plus avec la bibliothèque standard...

"""

# ╔═╡ 26930234-689a-4766-b33c-8260ee01e29f
cg(Diagonal(dA), db)

# ╔═╡ Cell order:
# ╟─6dec645a-8a0b-49c8-8ae9-37c7c7f2b25e
# ╠═c0340c58-c233-11eb-32fd-ddb8506d9ecb
# ╟─157ab7ce-ed31-430e-a230-438f02630753
# ╠═da7cad0e-40eb-404f-afdf-f3e1f5aa2307
# ╠═23e036d6-c67d-4e11-9d17-b87f6f0657c8
# ╟─5029cc06-922e-4c74-953d-cc660bb551c2
# ╠═9fd097d5-1696-4426-a2a0-e78e7d611df5
# ╟─a5cf0863-b6a7-4dab-ae80-68f58484bf17
# ╠═90f90d71-452f-45e0-bab9-3443b8e9cd98
# ╟─efac4950-d829-45b2-a9dc-958fb9d57996
# ╠═d46458b9-0b34-4fb7-bf94-442db95b9de1
# ╠═159d4b51-50e3-4955-833a-1510b321ec61
# ╟─e5ee30dd-633a-4c37-bf78-92a1a21d65cd
# ╠═c810af14-b553-4a25-9423-c18850dae690
# ╠═dbb335aa-b9b2-40e8-a431-a05a5639ea2b
# ╠═a61234f3-72d1-4d7a-b870-ef89542e7567
# ╠═dc49be41-b4b5-4fee-8e7f-c45969728e71
# ╠═eeaf4a3e-51f5-4951-a0a5-0f1e6c16c2cd
# ╟─38c3e98c-4cbc-4559-982c-6b9f94574ad2
# ╟─b1b12b1b-1cc4-4c1d-9352-95551d1b0aa3
# ╠═626f475f-bc14-4f82-bac8-3976bf3f25d0
# ╠═aa75ad03-74d9-4602-9219-e641f28a3a15
# ╠═64b6cc94-3297-4a89-8c6f-23f9da1e110d
# ╠═cb3eb1f1-c1e3-4b08-8243-e8d523cc355c
# ╠═b8f68ab8-84eb-4a11-bb5c-35d8b5f95057
# ╟─f61420f4-024e-4bcf-99dc-881ccccc5404
# ╠═42bbbcd8-1d50-42df-a700-42743d10e417
# ╠═26930234-689a-4766-b33c-8260ee01e29f
