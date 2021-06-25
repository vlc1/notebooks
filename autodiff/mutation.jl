### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 748a0783-bb0e-497c-ab10-3b1b03a02112
using ChainRulesCore, ChainRules

# ╔═╡ 65854f7a-25a7-484c-a1ef-fe459684c639
md"""
Re-read the [documentation](https://juliadiff.org/ChainRulesCore.jl/dev/index.html).

# Immutable types

"""

# ╔═╡ d32f7968-6cab-44d6-b8d0-83d1cf8e4b66
f(x) = x ^ 2

# ╔═╡ 87b62049-9860-4137-b644-6af0e9219e03
function ChainRulesCore.rrule(::typeof(f), x)
	y = f(x)
	function fpullback(Δy)
		return NO_FIELDS, 2x * Δy
	end
	y, fpullback
end

# ╔═╡ 920c51cb-888e-46e2-bbc8-8c603a07c554
function ∇f(x)
	y, pullback = rrule(f, x)
	_, y̅ = pullback(one(x))
	y̅
end

# ╔═╡ 31891d44-fedf-49be-80a1-445a164cf60b
∇f(2)

# ╔═╡ 4e7768f7-f317-434f-af22-ecb94ae6ec54
md"""
# Mutable types

"""

# ╔═╡ 481bee55-b07e-4096-a9d3-c3bf79c2b7b7
md"""
Probably the simplest type of mutable object is the 0-dimensional array.

"""

# ╔═╡ 00bd99f9-a9f1-4169-bc9a-98e6e9f77391
f!(y, x) = @. y = f(x)

# ╔═╡ c945b33a-0b4d-4ccf-820f-e923a028a602
function foo()
	x = rand(Float64, ())
	y = similar(x)
	f!(y, x)
	x, y
end

# ╔═╡ 53605702-c2de-11eb-05b3-9f69cde3f559
foo()

# ╔═╡ a55da757-15d0-4321-b376-320d486ebeef
md"""
Next step is to define the `rrule` for `f!`.

"""

# ╔═╡ Cell order:
# ╠═748a0783-bb0e-497c-ab10-3b1b03a02112
# ╟─65854f7a-25a7-484c-a1ef-fe459684c639
# ╠═d32f7968-6cab-44d6-b8d0-83d1cf8e4b66
# ╠═87b62049-9860-4137-b644-6af0e9219e03
# ╠═920c51cb-888e-46e2-bbc8-8c603a07c554
# ╠═31891d44-fedf-49be-80a1-445a164cf60b
# ╟─4e7768f7-f317-434f-af22-ecb94ae6ec54
# ╟─481bee55-b07e-4096-a9d3-c3bf79c2b7b7
# ╠═00bd99f9-a9f1-4169-bc9a-98e6e9f77391
# ╠═c945b33a-0b4d-4ccf-820f-e923a028a602
# ╠═53605702-c2de-11eb-05b3-9f69cde3f559
# ╟─a55da757-15d0-4321-b376-320d486ebeef
