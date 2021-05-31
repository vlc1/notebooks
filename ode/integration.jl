### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ c42a3e90-8250-11eb-2e52-174b1b7a87bc
using Plots

# ╔═╡ 3658aed6-824f-11eb-049c-cb2ce9117e65
md"""
Let us consider the solution of the following form of initial value problem
```math
\left \{ \begin{aligned}
\dot{y} \left ( t \right ) & = f \left [ t, y \left ( t \right ) \right ], \quad \forall t > 0, \\
y \left ( 0 \right ) & = y_0.
\end{aligned} \right .
```

As an example, let us consider the case where
```math
\left \{ \begin{aligned}
f & \colon \left ( t , y \right ) \mapsto t + 2 - y, \\
y_0 & = 0,
\end{aligned} \right .
```
the exact solution then reads
```math
y \colon t \mapsto t + 1 - \exp \left (-t \right ).
```

Let us know integrate the ODE over an interval ``I_n \equiv \left ( t_n, t_{n + 1} \right ]``. Stokes' theorem yields
```math
y \left ( t_{n + 1} \right ) - y \left ( t_n \right ) = \int_{t_n}^{t_{n + 1}} f \left [ t, y \left ( t \right ) \right ] \mathrm{d} t.
```

Let us know visualize the graph of the integrand in the right-hand side, namely
```math
g \colon t \mapsto f \left [ t, y \left ( t \right ) \right ].
```

"""

# ╔═╡ c6add1a4-8250-11eb-393d-0fe1c50dafef
y(t) = t + 1 - exp(-t)

# ╔═╡ d585e066-8250-11eb-242a-3fc27732c4ac
f(t, y) = t + 2 - y

# ╔═╡ e8638898-8250-11eb-1e62-fb168f780ad2
g(t) = f(t, y(t))

# ╔═╡ b6b62c80-830c-11eb-1707-b983785f3cb7
default(lw = 2)

# ╔═╡ f307d86c-8250-11eb-01a6-2fe25dc9e9bb
begin
	xlim, ylim = (0, 2), (-0.1, 2.1)
	x = range(xlim..., step = 0.01)

	fig = plot(; xlim, ylim)
	plot!(fig, x, g.(x), fillrange = zero.(x), fillalpha = 0.35, label = "g")
end

# ╔═╡ 1b786106-8254-11eb-31db-fb4ef030eb1a
md"""
# Forward Euler

"""

# ╔═╡ 145aa8ec-8255-11eb-3e5d-23a0e9a3b983
begin
	local x = 0.0:0.1:0.1
	plot!(fig, x, fill(g(0.0), size(x)), fillrange = zero.(x), fillalpha = 0.35, label = "")
end

# ╔═╡ Cell order:
# ╟─3658aed6-824f-11eb-049c-cb2ce9117e65
# ╠═c6add1a4-8250-11eb-393d-0fe1c50dafef
# ╠═d585e066-8250-11eb-242a-3fc27732c4ac
# ╠═e8638898-8250-11eb-1e62-fb168f780ad2
# ╠═c42a3e90-8250-11eb-2e52-174b1b7a87bc
# ╠═b6b62c80-830c-11eb-1707-b983785f3cb7
# ╠═f307d86c-8250-11eb-01a6-2fe25dc9e9bb
# ╟─1b786106-8254-11eb-31db-fb4ef030eb1a
# ╠═145aa8ec-8255-11eb-3e5d-23a0e9a3b983
