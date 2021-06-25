### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 38af4e9c-67f1-11eb-0d06-39beab11e515
using ChainRulesCore

# ╔═╡ c82aac64-672f-11eb-2e99-abd7cb5bbb60
md"""
From `Zygote`'s [documentation] : the gradient of a `getindex`
```julia
y = x[i...]
```
is a `setindex!`
```julia
x̄[i...] = ȳ.
```

`hook`, `checkpoint`...

Could be nice places to start :

1. <https://arxiv.org/abs/1810.07951>,
1. <https://arxiv.org/abs/1907.07587>.

Also check

1. `Zygote`'s internals (Mike Innes' notebooks),
1. If mutation/arrays too difficult, start with checkpointing (see `checkpointed`).

"""

# ╔═╡ 1136c1d4-67f0-11eb-3cd4-4d33fa845067
begin
	mutable struct Foo{T}
		val::T
	end

	Foo(::Type{T}) where {T} = Foo(rand(T))
end

# ╔═╡ 4c321e6e-67f0-11eb-1939-a73c4ca55cad
Base.getindex(x::Foo) = x.val

# ╔═╡ 6099fd86-67f0-11eb-12e3-39dfc087e5b8
Base.setindex!(x::Foo, val) = x.val = val

# ╔═╡ 724b0d2c-67f0-11eb-1f38-e9906e48705d
begin
	foo = Foo(1.0)
	foo[] = 2
	foo
end

# ╔═╡ 73401c02-67f4-11eb-02cb-737660c26398
Foo(Int)

# ╔═╡ 27b99738-67f0-11eb-0c0d-5331a43a7383
md"""
```math
f \colon \left ( x_1, x_2 \right ) \mapsto \mathtt{Foo(}x_1\mathtt{)[]} ^ 2 + x_2.
```

"""

# ╔═╡ fdbbbe38-67f0-11eb-0ec9-0372ced9ce00
function f(x₁, x₂)
	x₃ = Foo(x₁)
	x₄ = x₃[]
	x₅ = x₄ ^ 2
	x₆ = x₂ + x₅
end

# ╔═╡ 17ca8994-67f1-11eb-22aa-67c2e1e24990
x⃗ = x₁, x₂ = 4, 9

# ╔═╡ 47df0e84-67f1-11eb-1ba8-b57e8a3cf262
@scalar_rule(x + y, (One(), One()))

# ╔═╡ 62b56122-67f1-11eb-1744-c75546979298
@scalar_rule(x ^ y, (y * Ω / x, Ω * log(x)))

# ╔═╡ f3cd8482-67f1-11eb-256c-3d5b45caaaf8
@scalar_rule(x * y, (y, x))

# ╔═╡ Cell order:
# ╟─c82aac64-672f-11eb-2e99-abd7cb5bbb60
# ╠═1136c1d4-67f0-11eb-3cd4-4d33fa845067
# ╠═4c321e6e-67f0-11eb-1939-a73c4ca55cad
# ╠═6099fd86-67f0-11eb-12e3-39dfc087e5b8
# ╠═724b0d2c-67f0-11eb-1f38-e9906e48705d
# ╠═73401c02-67f4-11eb-02cb-737660c26398
# ╟─27b99738-67f0-11eb-0c0d-5331a43a7383
# ╠═fdbbbe38-67f0-11eb-0ec9-0372ced9ce00
# ╠═17ca8994-67f1-11eb-22aa-67c2e1e24990
# ╠═38af4e9c-67f1-11eb-0d06-39beab11e515
# ╠═47df0e84-67f1-11eb-1ba8-b57e8a3cf262
# ╠═62b56122-67f1-11eb-1744-c75546979298
# ╠═f3cd8482-67f1-11eb-256c-3d5b45caaaf8
