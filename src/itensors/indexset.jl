
"""
    hasinds(is)

Returns an anonymous function `x -> hasinds(x, is)` which
accepts an ITensor or IndexSet and returns `true` if the
ITensor or IndexSet has the indices `is`.
"""
ITensors.hasinds(s) = x -> hasinds(x, s)

@eval struct Order{x}
  (OrderT::Type{ <: Order})() = $(Expr(:new, :OrderT))
end

Order(x) = Order{x}()

# TODO: make IndexSet{N} type stable for element type
# and storage type.
Base.filter(O::Order{N},
            f::Function,
            is::IndexSet) where {N} = IndexSet{N}(filter(f,
                                                         Tuple(is)))

Base.filter(O::Order,
            is::IndexSet,
            args...; kwargs...) = filter(O,
                                         ITensors.fmatch(args...;
                                                         kwargs...),
                                         is)

# TODO: export this. Maybe make this more
# specific than a typedef?
const filterinds = inds

# intersect
ITensors.commoninds(::Order{N},
           A...;
           kwargs...) where {N} =
  IndexSet{N}(intersect(ITensors.itensor2inds.(A)...;
                        kwargs...)...)

# symdiff
ITensors.noncommoninds(::Order{N},
              A...;
              kwargs...) where {N} =
  IndexSet{N}(symdiff(ITensors.itensor2inds.(A)...;
                      kwargs...)...)

# setdiff
ITensors.uniqueinds(::Order{N},
           A...;
           kwargs...) where {N} =
  IndexSet{N}(setdiff(ITensors.itensor2inds.(A)...;
                      kwargs...)...)

# union
ITensors.unioninds(::Order{N},
          A...;
          kwargs...) where {N} =
  IndexSet{N}(union(ITensors.itensor2inds.(A)...;
                    kwargs...)...)

