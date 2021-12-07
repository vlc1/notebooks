### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ d6b4b4d8-53af-11ec-2309-f383ef2c77dc
using Enzyme

# ╔═╡ a63502f2-317b-4d07-8009-6b9639907bf5
md"""
JuliaCon 2021 [talk](https://www.youtube.com/watch?v=gFfePK44ICk) on `Enzyme.jl`.

"""

# ╔═╡ ea8aee31-bd41-491e-8b0f-8cf9dc07dd75
square(x) = x * x

# ╔═╡ f3c8468f-474f-419d-a411-08709d114798
foo(f, x) = x .* f.(x)

# ╔═╡ 147992a0-7f97-4653-8611-b2692c0d0f2b
md"""
# Activity of constant arguments (`Const`)

Default activity for values is `Const`.

"""

# ╔═╡ 2ee3f10e-505b-4c8b-ad1f-c4c99a676f20
autodiff(square, 1.0)

# ╔═╡ 846652f9-9fb4-447d-a26b-6aedd9eabc66
autodiff(foo, sqrt, Active(2))

# ╔═╡ f39f960f-ba16-4001-8e9b-8a8b3191fa57
md"""
# Activity of scalar/immutable arguments

"""

# ╔═╡ 99602d71-6eac-4fa4-a3a6-9da764ac7134
autodiff(square, Active(1.0))

# ╔═╡ 9d8cae6d-0815-4fab-9839-25cbfdcb2691
md"""
# Activity of mutable arguments (`Duplicated`/`DuplicatedNoNeed`)

"""

# ╔═╡ ee5ae350-5481-4d97-8a95-447bcb87aaa9
function square!(y, x)
	y[] = x[] ^ 2
	nothing
end

# ╔═╡ e9bbae72-13fa-4229-a2f9-5f132614f13f
let
	x = Ref(4.0)
	dx = Ref(0.0) # accumulates

	y = Ref(-4.0) # overwritten when Duplicated
	dy = Ref(1.0) # passed onto

	autodiff(square!, Duplicated(y, dy), Duplicated(x, dx))

	y[], dy[], x[], dx[]
end

# ╔═╡ 3e57bce7-2462-49fc-92d3-4286b1fa7a6c
function cube!(y, x)
	y .= x .^ 2
	nothing
end

# ╔═╡ ea3bdead-fa30-42e6-9d80-5dbe3a3202f0
let
	x = rand(4, 8)
	dx = fill!(similar(x), 0)

	y = rand(4, 8)
	dy = fill!(similar(y), 1)

	autodiff(cube!, Duplicated(y, dy), Duplicated(x, dx))

	x, dx, y, dy
end

# ╔═╡ 88584215-189a-4bb9-9ed1-b1b2cfa301bd
md"""
How to compose?

"""

# ╔═╡ 9a5ab937-5604-4c93-b06a-cc05f015a688
h = √

# ╔═╡ 23c43e4a-9982-4465-91e3-0a805033e42c
function g!(z, x)
	z .= h.(x)
	nothing
end

# ╔═╡ f2d42244-0dd1-46fe-a720-00a8ab9e857e
function f!(y, x)
	z = similar(x)
	g!(z, x)
	y .= z
	nothing
end

# ╔═╡ da9a262b-79c6-4de6-b65b-c701bd5ce466
let
	x = rand(4, 6)
	dx = zero(x) # accumulates

	y = Array(x) # overwritten when Duplicated
	dy = fill!(similar(y), 1) # passed onto

	autodiff(f!, Duplicated(y, dy), Duplicated(x, dx))
end

# ╔═╡ fb0cccf0-2278-41d2-883a-e110e443fd11
x = rand(4, 5)

# ╔═╡ fa762cf7-5b07-4c46-9cfa-9f3971bb8d65
fill!(x, 1)

# ╔═╡ 1d76741c-9552-4d22-a5ae-1df69a6d2485
fill(x, 1)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[compat]
Enzyme = "~0.7.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-rc3"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Enzyme]]
deps = ["Adapt", "CEnum", "Enzyme_jll", "GPUCompiler", "LLVM", "Libdl", "ObjectFile"]
git-tree-sha1 = "b1beaca0a3d99bf62fb9ea18d287d3ee29dc9782"
uuid = "7da242da-08ed-463a-9acd-ee780be4f1d9"
version = "0.7.2"

[[deps.Enzyme_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0c1182e39da5e26fec00a4af9a63ab59fc175d8e"
uuid = "7cc45869-7501-5eee-bdea-0790c847d4ef"
version = "0.0.23+0"

[[deps.ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[deps.GPUCompiler]]
deps = ["ExprTools", "InteractiveUtils", "LLVM", "Libdl", "Logging", "TimerOutputs", "UUIDs"]
git-tree-sha1 = "6cf994358b3821ea446c43dea08c38aceb60a0cc"
uuid = "61eb1bfa-7361-4325-ad38-22787b887f55"
version = "0.13.8"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "7cc22e69995e2329cc047a879395b2b74647ab5f"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.7.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c5fc4bef251ecd37685bea1c4068a9cfa41e8b9a"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.13+0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.ObjectFile]]
deps = ["Reexport", "StructIO"]
git-tree-sha1 = "55ce61d43409b1fb0279d1781bf3b0f22c83ab3b"
uuid = "d8793406-e978-5875-9003-1fc021f44a92"
version = "0.3.7"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.StructIO]]
deps = ["Test"]
git-tree-sha1 = "010dc73c7146869c042b49adcdb6bf528c12e859"
uuid = "53d494c1-5632-5724-8f4c-31dff12d585f"
version = "0.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "7cb456f358e8f9d102a8b25e8dfedf58fa5689bc"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.13"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─a63502f2-317b-4d07-8009-6b9639907bf5
# ╠═d6b4b4d8-53af-11ec-2309-f383ef2c77dc
# ╠═ea8aee31-bd41-491e-8b0f-8cf9dc07dd75
# ╠═f3c8468f-474f-419d-a411-08709d114798
# ╟─147992a0-7f97-4653-8611-b2692c0d0f2b
# ╠═2ee3f10e-505b-4c8b-ad1f-c4c99a676f20
# ╠═846652f9-9fb4-447d-a26b-6aedd9eabc66
# ╟─f39f960f-ba16-4001-8e9b-8a8b3191fa57
# ╠═99602d71-6eac-4fa4-a3a6-9da764ac7134
# ╟─9d8cae6d-0815-4fab-9839-25cbfdcb2691
# ╠═ee5ae350-5481-4d97-8a95-447bcb87aaa9
# ╠═e9bbae72-13fa-4229-a2f9-5f132614f13f
# ╠═3e57bce7-2462-49fc-92d3-4286b1fa7a6c
# ╠═ea3bdead-fa30-42e6-9d80-5dbe3a3202f0
# ╟─88584215-189a-4bb9-9ed1-b1b2cfa301bd
# ╠═9a5ab937-5604-4c93-b06a-cc05f015a688
# ╠═23c43e4a-9982-4465-91e3-0a805033e42c
# ╠═f2d42244-0dd1-46fe-a720-00a8ab9e857e
# ╠═da9a262b-79c6-4de6-b65b-c701bd5ce466
# ╠═fb0cccf0-2278-41d2-883a-e110e443fd11
# ╠═fa762cf7-5b07-4c46-9cfa-9f3971bb8d65
# ╠═1d76741c-9552-4d22-a5ae-1df69a6d2485
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
