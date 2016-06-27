module RecursiveFiltering

using KernelDensityEstimate, TransformUtils

export
  SysState,
  filterpropagate!,
  filterupdate!

include("KalmanFilter.jl")
include("KDEFilter.jl")

end
