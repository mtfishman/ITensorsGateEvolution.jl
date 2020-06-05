#
# Op
#

struct Op{N, Params <: NamedTuple}
  name::String
  sites::NTuple{N,Int}
  params::Params
end

Op(name::String, sites::Tuple{Vararg{Int}}) =
  Op(name, sites, NamedTuple())

Op(name::String, site::Int,
   params::NamedTuple = NamedTuple()) =
  Op(name, (site,), params)

Op(name::String, sites::Int...) =
  Op(name, sites)

Base.convert(::Type{Op}, t::Tuple) = Op(t...)

function ITensors.op(sites::Vector{<:Index},
                     o::Op)
  return op(sites, o.name, o.sites...; o.params...)
end

function Base.show(io::IO, o::Op)
  print(io, "\"$(o.name)\"")
  if length(o.sites) == 1
    print(io, "($(only(o.sites)))")
  else
    print(io, o.sites)
  end
  if !isempty(o.params)
    print(io, " ")
    if length(o.params) == 1
      print(io, "($(only(keys(o.params))) = $(only(o.params)))")
    else
      print(io, o.params)
    end
  end
end

#
# ProductOps
#

struct ProductOps
  data::Vector{Op}
end

ProductOps() = ProductOps(Op[])

Base.copy(C::ProductOps) = ProductOps(copy(C.data))

Base.push!(C::ProductOps, O) =
  (push!(C.data, O);
   return C)

Base.pushfirst!(C::ProductOps, O) =
  (pushfirst!(C.data, O);
   return C)

Base.iterate(C::ProductOps, args...) = iterate(C.data, args...)
Base.length(C::ProductOps) = length(C.data)

import Base: *

C::ProductOps * O = push!(copy(C), O)

Op * C::ProductOps = pushfirst!(copy(C), O)

function Base.show(io::IO, C::ProductOps)
  println("ProductOps")
  for o in C.data
    println(io, o)
  end
end

#
# qubit
#

ITensors.siteinds(::TagType"qubit", N::Int) =
  [Index(2, "Site,qubit,n=$n") for n in 1:N]

function ITensors.state(::TagType"qubit",
                        st::AbstractString)
  if st == "0"
    return 1
  elseif st == "1"
    return 2
  end
  throw(ArgumentError("State string \"$st\" not recognized for SpinHalf site"))
  return 0
end

function ITensors.op(::TagType"qubit",
                     s::Index,
                     opname::AbstractString;
                     kwargs...)::ITensor{2}
  Up = s(1)
  UpP = s'(1)
  Dn = s(2)
  DnP = s'(2)

  Op = emptyITensor(s',dag(s))

  if opname == "I"
    Op[UpP, Up] = 1
    Op[DnP, Dn] = 1
  elseif opname == "X"
    Op[UpP, Dn] = 1
    Op[DnP, Up] = 1
  elseif opname == "iY"
    Op[UpP, Dn] = 1
    Op[DnP, Up] = -1
  elseif opname == "Y"
    Op = complex(Op)
    Op[UpP, Dn] = -im
    Op[DnP, Up] = im
  elseif opname == "Z"
    Op[UpP, Up] = 1
    Op[DnP, Dn] = -1
  elseif opname == "H"
    Op[UpP, Up] = 1/sqrt(2)
    Op[UpP, Dn] = 1/sqrt(2)
    Op[DnP, Up] = 1/sqrt(2)
    Op[DnP, Dn] = -1/sqrt(2)
  else
    throw(ArgumentError("Operator name '$opname' not recognized for SpinHalfSite"))
  end
  return Op
end

#
# product
#

"""
    ops(::Vector{Index}, C::ProductOps)

Return a Vector of ITensors corresponding to the input circuit.
"""
function ops(s::Vector{<:Index}, C::ProductOps)
  return [op(s, c) for c in C]
end

"""
    movesite(::MPS, n1::Int, n2::Int)

Create a new MPS where the site at `n1` is moved to `n2`
"""
function movesite(ψ::MPS, n1::Int, n2::Int)
  r = n1:n2-1
  ortho = "left"
  if n1 > n2
    r = reverse(n2:n1-1)
    ortho = "right"
  end
  for n in r
    ψ = swapbondsites(ψ, n; ortho = ortho)
  end
  return ψ
end

"""
    sitedict(::Vector{<:Index})

Return a dictionary that maps a Vector of indices to
the integer position of the Index in the Vector.
"""
function sitedict(sites::Vector{IndexT}) where {IndexT <: Index}
  d = Dict{IndexT, Int}()
  for (n, s) in enumerate(sites)
    d[s] = n
  end
  return d
end

"""
    product(ops::Vector{<:ITensor}, ψ::MPS)

Apply the ITensors `ops` to the MPS `ψ`.
"""
function ITensors.product(ops::Vector{<:ITensor},
                          ψ0::MPS)
  ψ = copy(ψ0)
  d = sitedict(siteinds(ψ0))
  for o in ops
    n = d[firstind(o; plev = 0)]
    orthogonalize!(ψ, n)
    ψ[n] = noprime(o * ψ[n])
  end
  return ψ
end

const apply = product

