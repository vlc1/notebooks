### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ ecf08d34-777d-4743-8189-ce570b3a89e2
using Plots

# ╔═╡ ecf3120b-fc6f-4884-9f19-dd3e8e6ecdfa
using BenchmarkTools

# ╔═╡ 3da4102e-acdc-11eb-08c8-dbae38bb2fea
md"""
# Divided differences

"""

# ╔═╡ 121ff2ec-6a36-4c17-80a0-749a6c634e1a
function mdd!(ys, xs)
	n = length(ys)

	n == 1 && return

	m = length(xs)
	j = m - n

	for i in reverse(2:n)
		ys[i] = (ys[i] - ys[i - 1]) / (xs[i + j] - xs[i - 1])
	end

	mdd!(view(ys, 2:n), xs)
end

# ╔═╡ a96d94f8-59f7-495f-9c5c-24db02656067
function hoener(zs, xs, x)
	n = length(zs)

	n == 1 && return last(zs)

	first(zs) + (x - first(xs)) * hoener(view(zs, 2:n), view(xs, 2:n), x)
end

# ╔═╡ 7b6e8053-1d49-4127-bd9a-1779368291ef
function hoener2(zs, xs, x)
	n, y = length(zs), last(zs)
	for i in reverse(1:n - 1)
		y = zs[i] + (x - xs[i]) * y
	end
	y
end

# ╔═╡ 4235b72b-5e90-4f85-a4d4-dc4a9cde35c6
hoener(zs, xs) = interpolant(x) = hoener(zs, xs, x)

# ╔═╡ be9bbacf-0d92-4a67-89bb-9cb167e05312
md"""
## Application

"""

# ╔═╡ 1be18d7f-7a60-40a7-b8e0-19b3a24309b4
f(x, α = 1 / √8) = α ^ 2 / (x ^ 2 + α ^ 2)

# ╔═╡ 6180c802-c04a-4158-99fb-6c78426b7522
md"""
### Uniform spacing

"""

# ╔═╡ 42b1bd32-7182-4a11-b927-38facc97b09c
begin
	local n = 6
	xs = -1:2 / n:1
	ys = f.(xs)
	zs = copy(ys)
	mdd!(zs, xs)
end

# ╔═╡ 091c4a61-b6c8-47b6-a279-2b07758fd5dc
begin
	local fig = plot(
		xlim = (-1.1, 1.1),
		ylim = (-0.1, 1.1)
	)
	scatter!(fig, xs, ys, label = "data")
	plot!(fig, f, label = "exact")
	plot!(fig, hoener(zs, xs), label = "uniform")
end

# ╔═╡ f672ff8b-ac4a-43fd-b48a-ab8669daa940
md"""
### Tchebychev points

"""

# ╔═╡ 66cf2413-a46a-4ecc-ac18-e9d9c4daabcb
begin
	local n = 6
	xt = [cos((2i - 1)π / (2n + 2)) for i in 1:n + 1]
	yt = f.(xt)
	zt = copy(yt)
	mdd!(zt, xt)
end

# ╔═╡ 4aef7add-47ca-440e-926a-30a3d695ab14
begin
	local fig = plot(
		xlim = (-1.1, 1.1),
		ylim = (-0.1, 1.1)
	)
	scatter!(fig, xt, yt, label = "data")
	plot!(fig, f, label = "exact")
	plot!(fig, hoener(zt, xt), label = "Tchebychev")
end

# ╔═╡ 866c9be0-a164-403e-8751-ea88fa6b2d00
@btime hoener2($zs, $xs, 0.)

# ╔═╡ Cell order:
# ╟─3da4102e-acdc-11eb-08c8-dbae38bb2fea
# ╠═121ff2ec-6a36-4c17-80a0-749a6c634e1a
# ╠═a96d94f8-59f7-495f-9c5c-24db02656067
# ╠═7b6e8053-1d49-4127-bd9a-1779368291ef
# ╠═4235b72b-5e90-4f85-a4d4-dc4a9cde35c6
# ╟─be9bbacf-0d92-4a67-89bb-9cb167e05312
# ╠═1be18d7f-7a60-40a7-b8e0-19b3a24309b4
# ╟─6180c802-c04a-4158-99fb-6c78426b7522
# ╠═42b1bd32-7182-4a11-b927-38facc97b09c
# ╠═ecf08d34-777d-4743-8189-ce570b3a89e2
# ╠═091c4a61-b6c8-47b6-a279-2b07758fd5dc
# ╟─f672ff8b-ac4a-43fd-b48a-ab8669daa940
# ╠═66cf2413-a46a-4ecc-ac18-e9d9c4daabcb
# ╠═4aef7add-47ca-440e-926a-30a3d695ab14
# ╠═ecf3120b-fc6f-4884-9f19-dd3e8e6ecdfa
# ╠═866c9be0-a164-403e-8751-ea88fa6b2d00
