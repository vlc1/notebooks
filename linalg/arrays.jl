struct MyArray{T,N,AA<:AbstractArray{T,N}} <: AbstractArray{T,N}
    parent::AA
end

import Base: parent, size, getindex, setindex!, similar

parent(A::MyArray) = getproperty(A, :parent)

size(A::MyArray) = size(parent(A))

getindex(A::MyArray, i...) = getindex(parent(A), i...)

setindex!(A::MyArray, v, i...) = setindex!(parent(A), v, i...)

similar(A::MyArray, ::Type{T}, dims::Dims) where {T} =
    MyArray(similar(parent(A), T, dims))

