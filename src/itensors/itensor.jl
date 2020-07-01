
Base.eachindex(A::ITensor) = eachindex(ITensors.tensor(A))

Base.setindex!(T::ITensor, x::Number, I::CartesianIndex) =
  setindex!(T, x, Tuple(I)...)

#TODO: Make a constructor
#  itensor(::Vector{Pair{Block, Array}}, ::IndexSet)
#for example:
#  itensor([Block(1,1) => randn(2,2),
#           Block(2,2) => randn(3,3)], i', dag(i)))
#to set nonzero blocks.
function ITensors.itensor(A::Array{ElT},
                          inds::ITensors.QNIndexSet) where {ElT <: Number}
  length(A) ≠ dim(inds) && throw(DimensionMismatch("In ITensor(::Array, ::IndexSet), length of Array ($(length(A))) must match total dimension of IndexSet ($(dim(inds)))"))
  T = emptyITensor(ElT, inds)
  A = reshape(A, dims(inds))
  for vs in eachindex(T)
    Avs = A[vs]
    if !iszero(Avs)
      T[vs] = A[vs]
    end
  end
  return T
end

ITensors.itensor(A::Array{<:Number},
                 inds::ITensors.QNIndex...) =
  itensor(A, IndexSet(inds...))

# TODO: make this function
#IndexSet{N,Index{Array{Pair{QN,Int64},1}},Tuple{Vararg{Index{Array{Pair{QN,Int64},1}},N}}} where N(::Index{Array{Pair{QN,Int64},1}}, ::Index{Array{Pair{QN,Int64},1}})

function hascommoninds(A::ITensor,
                       B::ITensor; kwargs...)
  return !isnothing(commonind(A, B; kwargs...))
end

