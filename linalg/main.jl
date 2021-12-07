using LinearAlgebra

n = 4

# Julia's native arrays
A = Tridiagonal(-ones(n-1), 2ones(n), -ones(n-1))
b = rand(n)
x = similar(b)

# Wrapping with our own type
include("arrays.jl")

A̅ = MyArray(A)
b̅ = MyArray(b)
x̅ = similar(b̅)

using IterativeSolvers

cg!(x, A, b)
cg!(x̅, A̅, b̅)

x, x̅

