module RecursiveFiltering

using KernelDensityEstimate

export
  SysState,
  filterpropagate!,
  filterupdate!

include("KalmanFilter.jl")
include("KDEFilter.jl")

end
