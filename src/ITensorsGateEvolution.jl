module ITensorsGateEvolution

using ITensors

export ProductOps,
       ops,
       apply

# Extensions to ITensors
include("itensors/tupletools.jl")
include("itensors/index.jl")
include("itensors/indexset.jl")
include("itensors/itensor.jl")
include("itensors/mps.jl")

# ProductOps type
include("productops.jl")

# Qubit site type
include("qubit.jl")

# Apply function
include("apply.jl")

end
