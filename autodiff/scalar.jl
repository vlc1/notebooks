### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 1bdb5f04-6666-11eb-39e9-e30a88d1d554
using ChainRulesCore

# ╔═╡ 40f2b83a-6726-11eb-311e-25fd12a7e3d0
using BenchmarkTools, Zygote

# ╔═╡ 2b1ba488-6666-11eb-3227-15af872ab470
md"""
Let's get our own AD implemention using `ChainRulesCore.jl` working on Nocedal and Wright's example :
```math
f \colon \left ( x_1, x_2, x_3 \right ) \mapsto \frac{x_1 x_2 \sin \left ( x_3 \right ) + \exp \left ( x_1 x_2 \right )}{x_3}.
```

To do that we'll need to specify a few basic rules (check `ChainRules.jl`'s `src/rulesets` directory).

"""

# ╔═╡ aa767f40-6667-11eb-3cb0-3b2fb0dbdadf
x⃗ = (x₁, x₂, x₃) = 1.0, 2.0, π / 2

# ╔═╡ f716ab12-6723-11eb-2c26-8d9c007bb602
y = (4 + 2exp(2)) / π

# ╔═╡ 812f9fdc-6666-11eb-31fa-eb63ec4d5abb
begin
	@scalar_rule(sin(x), cos(x)) # more efficient implementation in ChainRules
	@scalar_rule(exp(x), Ω)
	@scalar_rule(x + y, (One(), One()))
	@scalar_rule(x * y, (one(x) * y, x * one(y)))
	@scalar_rule(x / y, (one(x) / y, -(Ω / y)))
end

# ╔═╡ 77736092-6723-11eb-0b74-07ed20366007
g⃗ = 4(exp(2) + 1) / π, 2(exp(2) + 1) / π, -4(2 + exp(2)) / π ^ 2

# ╔═╡ 50a9444c-6668-11eb-36fd-0fa0f33468da
md"""
Nocedal and Wright evaluate ``f`` as follows (see graphe in their book) :

"""

# ╔═╡ 206b9b4e-671d-11eb-35f4-9da75bb02bc8
function f(x₁, x₂, x₃)
	x₄ = x₁ * x₂
	x₅ = sin(x₃)
	x₆ = exp(x₄)
	x₇ = x₄ * x₅
	x₈ = x₆ + x₇
	x₉ = x₈ / x₃
	x₉
end

# ╔═╡ 653c7d7e-671d-11eb-2959-a9d25462b767
md"""
And we can easily check that it works by comparing the numerical and analytical value.

"""

# ╔═╡ 3415e096-671d-11eb-2a9a-9be735949d7e
f(x⃗...), y

# ╔═╡ 6869576c-6667-11eb-2e87-99a9bb078d00
md"""
## Forward mode

"""

# ╔═╡ c4ae2e80-6723-11eb-2ee6-311bb21f5158
function forward(x₁, x₂, x₃)
	x₉, ẋ₁ = forward(x₁, x₂, x₃, one(x₁), zero(x₂), zero(x₃))
	x₉, ẋ₂ = forward(x₁, x₂, x₃, zero(x₁), one(x₂), zero(x₃))
	x₉, ẋ₃ = forward(x₁, x₂, x₃, zero(x₁), zero(x₂), one(x₃))

	x₉, (ẋ₁, ẋ₂, ẋ₃)
end

# ╔═╡ 65fbd406-6724-11eb-1179-319b36ed6cd0
function forward(x₁, x₂, x₃, ẋ₁, ẋ₂, ẋ₃)
	x₄, ẋ₄ = frule((Zero(), ẋ₁, ẋ₂), *, x₁, x₂)
	x₅, ẋ₅ = frule((Zero(), ẋ₃), sin, x₃)
	x₆, ẋ₆ = frule((Zero(), ẋ₄), exp, x₄)
	x₇, ẋ₇ = frule((Zero(), ẋ₄, ẋ₅), *, x₄, x₅)
	x₈, ẋ₈ = frule((Zero(), ẋ₆, ẋ₇), +, x₆, x₇)
	x₉, ẋ₉ = frule((Zero(), ẋ₈, ẋ₃), /, x₈, x₃)

	x₉, ẋ₉
end

# ╔═╡ 71cf923a-6667-11eb-247a-49df2199e459
md"""
## Reverse mode

"""

# ╔═╡ a15b9854-6668-11eb-2423-37a1471df011
function reverse(x₁, x₂, x₃)

	# forward sweep
	x₄, pullback₄ = rrule(*, x₁, x₂)
	x₅, pullback₅ = rrule(sin, x₃)
	x₆, pullback₆ = rrule(exp, x₄)
	x₇, pullback₇ = rrule(*, x₄, x₅)
	x₈, pullback₈ = rrule(+, x₆, x₇)
	x₉, pullback₉ = rrule(/, x₈, x₃)

	# reverse sweep
	x̅₉ = one(x₉)
	_, x̅₈, x̅₃ = pullback₉(x̅₉)
	_, x̅₆, x̅₇ = pullback₈(x̅₈)
	_, x̅₄, x̅₅ = pullback₇(x̅₇)
	_, ∂ = pullback₆(x̅₆); x̅₄ += ∂
	_, ∂ = pullback₅(x̅₅); x̅₃ += ∂
	_, x̅₁, x̅₂ = pullback₄(x̅₄)

	x₉, (x̅₁, x̅₂, x̅₃)
end

# ╔═╡ 421f47fa-6726-11eb-3e0c-75d9b9c195e3
md"""
# Wrapping it up

## Validation

Let's check first that the results are correct :

1. The function ``f`` itself,

"""

# ╔═╡ 3a45eaa2-672b-11eb-2fff-9bd486f6e6e6
f(x⃗...), y

# ╔═╡ 4b07fa24-672b-11eb-2898-6b2bd9ce32bb
md"""
2. And its gradient.

"""

# ╔═╡ fe824758-6724-11eb-18f9-2bfeb0781141
(y, g⃗), forward(x⃗...), reverse(x⃗...)

# ╔═╡ 7a432296-6726-11eb-2fe4-938fa0a8cdad
md"""
## Benchmark

Let's now benchmark the hard-coded implementation *vs.* `Zygote`'s.

"""

# ╔═╡ 55b4e8fc-672a-11eb-151e-73ca8b44bdbf
@benchmark f($x⃗...)

# ╔═╡ 58cb5818-6726-11eb-35c9-17744e8cdf5d
@benchmark forward($x⃗...)

# ╔═╡ 643e28ba-6726-11eb-38ff-33b3fd8c3f87
@benchmark reverse($x⃗...)

# ╔═╡ 94ea969e-6726-11eb-3ac7-dd24e5c077bc
@benchmark gradient(f, $x⃗...)

# ╔═╡ 4d07f134-672a-11eb-0212-19c7bfa9b36b
md"""
## Conclusions

1. Time elapsed for reverse-mode gradient evaluation is double the time required to evaluation ``f``.
1. `Zygote`'s gradient is as fast as the hard-coded reverse mode!
1. As expected (scalar-valued function), forward-mode gradient evaluation is slower (four times, here).

"""

# ╔═╡ Cell order:
# ╟─2b1ba488-6666-11eb-3227-15af872ab470
# ╠═1bdb5f04-6666-11eb-39e9-e30a88d1d554
# ╠═812f9fdc-6666-11eb-31fa-eb63ec4d5abb
# ╠═aa767f40-6667-11eb-3cb0-3b2fb0dbdadf
# ╠═f716ab12-6723-11eb-2c26-8d9c007bb602
# ╠═77736092-6723-11eb-0b74-07ed20366007
# ╟─50a9444c-6668-11eb-36fd-0fa0f33468da
# ╠═206b9b4e-671d-11eb-35f4-9da75bb02bc8
# ╟─653c7d7e-671d-11eb-2959-a9d25462b767
# ╠═3415e096-671d-11eb-2a9a-9be735949d7e
# ╟─6869576c-6667-11eb-2e87-99a9bb078d00
# ╠═c4ae2e80-6723-11eb-2ee6-311bb21f5158
# ╠═65fbd406-6724-11eb-1179-319b36ed6cd0
# ╟─71cf923a-6667-11eb-247a-49df2199e459
# ╠═a15b9854-6668-11eb-2423-37a1471df011
# ╟─421f47fa-6726-11eb-3e0c-75d9b9c195e3
# ╠═3a45eaa2-672b-11eb-2fff-9bd486f6e6e6
# ╟─4b07fa24-672b-11eb-2898-6b2bd9ce32bb
# ╠═fe824758-6724-11eb-18f9-2bfeb0781141
# ╟─7a432296-6726-11eb-2fe4-938fa0a8cdad
# ╠═40f2b83a-6726-11eb-311e-25fd12a7e3d0
# ╠═55b4e8fc-672a-11eb-151e-73ca8b44bdbf
# ╠═58cb5818-6726-11eb-35c9-17744e8cdf5d
# ╠═643e28ba-6726-11eb-38ff-33b3fd8c3f87
# ╠═94ea969e-6726-11eb-3ac7-dd24e5c077bc
# ╟─4d07f134-672a-11eb-0212-19c7bfa9b36b
