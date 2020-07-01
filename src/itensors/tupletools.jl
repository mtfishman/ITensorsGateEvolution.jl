
function Base.diff(t::NTuple{N}) where N
  return ntuple(i -> t[i+1] - t[i], Val(N-1))
end

