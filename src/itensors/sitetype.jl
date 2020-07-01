
# Version that takes `op("S=1/2", "X", Index(2, "i"))`.
# Shortcuts searching for the correct tag, so is a bit faster.
ITensors.op(st::AbstractString,
            on::AbstractString,
            s::Index...; kwargs...) = op(SiteType(st), OpName(on), s...; kwargs...)

