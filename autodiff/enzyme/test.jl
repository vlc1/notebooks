A = rand(10)
B = rand(10)
R = similar(A)

dA = zero(A)
dB = zero(B)
dR = fill!(similar(R), 1)

function foo!(R, A, B)
    R .= A .+ B
    nothing
end

autodiff(foo!, Const, Duplicated(R, dR), Duplicated(A, dA), Duplicated(B, dB))

