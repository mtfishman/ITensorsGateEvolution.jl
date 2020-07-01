#
# qubit
#

function ITensors.siteinds(::TagType"qubit", N::Int;
                           conserve_qns::Bool = false)
  if conserve_qns
    space = [QN() => 2]
  else
    space = 2
  end
  [Index(space, "Site,qubit,n=$n") for n in 1:N]
end

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

ITensors.op(::TagType"qubit",
            ::OpName"Id",
            s::Index) =
  itensor([1 0; 0 1], s', dag(s))

ITensors.op(::TagType"qubit",
            ::OpName"X",
            s::Index) = 
  itensor([0.0 1.0; 1.0 0.0], s', dag(s))

ITensors.op(::TagType"qubit",
            ::OpName"iY",
            s::Index) = 
  itensor([0 1; -1 0], s', dag(s))

ITensors.op(::TagType"qubit",
            ::OpName"Y",
            s::Index) = 
  itensor([0 -im; im 0], s', dag(s))

ITensors.op(::TagType"qubit",
            ::OpName"Z",
            s::Index) = 
  itensor([1 0; 0 -1], s', dag(s))

ITensors.op(::TagType"qubit",
            ::OpName"H",
            s::Index) = 
  itensor([1/sqrt(2) 1/sqrt(2);
           1/sqrt(2) -1/sqrt(2)], s', dag(s))

ITensors.op(::TagType"qubit",
            ::OpName"Rx",
            s::Index; θ::Number) =
  itensor([cos(θ/2)      -im*sin(θ/2.);
           -im*sin(θ/2.) cos(θ/2.)], s', dag(s))

ITensors.op(::TagType"qubit",
            ::OpName"Sw",
            s1::Index,
            s2::Index) =
  itensor([1 0 0 0;
           0 0 1 0;
           0 1 0 0;
           0 0 0 1],
          s1', s2',
          dag(s1), dag(s2))

ITensors.op(::TagType"qubit",
            ::OpName"rand",
            s::Index...) =
  randomITensor(prime.(s)...,
                dag.(s)...)

ITensors.op(::TagType"qubit",
            ::OpName"Cx",
            s1::Index,
            s2::Index) =
  itensor([1 0 0 0;
           0 1 0 0;
           0 0 0 1;
           0 0 1 0],
          s1', s2',
          dag(s1), dag(s2))

ITensors.op(::TagType"qubit",
            ::OpName"T",
            s1::Index,
            s2::Index,
            s3::Index) =
  itensor([1 0 0 0 0 0 0 0;
           0 1 0 0 0 0 0 0;
           0 0 1 0 0 0 0 0;
           0 0 0 1 0 0 0 0;
           0 0 0 0 1 0 0 0;
           0 0 0 0 0 1 0 0;
           0 0 0 0 0 0 0 1;
           0 0 0 0 0 0 1 0],
          s1', s2', s3',
          dag(s1), dag(s2), dag(s3))

ITensors.op(::TagType"qubit",
            ::OpName"noise",
            s::Index...;
            krausind = hasqns(s[1]) ? Index([QN() => 2], "kraus") : Index(2, "kraus")) =
  randomITensor(prime.(s)..., dag.(s)..., krausind)

