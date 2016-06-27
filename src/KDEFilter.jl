
# Work in progress here

# not going to work proper yet
function propagate(b1Xl::BallTreeDensity, bDxb1::Array{Float64,1}, s::Array{Float64,1})
  pts1 = getPoints(b1Xl)
  pts = zeros(size(pts1))
  for i in 1:size(pts,2)
    ent = [s[1]*randn();s[2]*randn();s[3]*randn()]
    b1Tb = inv(SE2(bDxb1+ent))
    pt = vec(pts1[:,i])
    bTl = SE2([pt[1];pt[2];0.0])
    pts[:,i] = vec(se2vee(b1Tb*bTl))[1:2]
  end
  return kde!( pts, "lcv");
end

# not going to work proper yet
function update(bhatXl::BallTreeDensity, z::Array{Float64,1}, s::Array{Float64,1}; N::Int=75)

  bXl = p2cPtsKDE(z,s, N=N)
  # take the product between predicted and measured position
  dummy = kde!(rand(2,N),[1.0])
  pGM, = prodAppxMSGibbsS(dummy, [bhatXl, bXl], Union{}, Union{}, 5)

  # error("update -- Not ready")
  return kde!(pGM, "lcv")
end
function updatelin(bhatXl::BallTreeDensity, z::Array{Float64,1}, s::Array{Float64,1}; N::Int=75)

  bXl = resample(kde!((z')',s),N)
  # take the product between predicted and measured position
  dummy = kde!(rand(length(z),N),[1.0])
  pGM, = prodAppxMSGibbsS(dummy, [bhatXl, bXl], Union{}, Union{}, 5)

  # error("update -- Not ready")
  return kde!(pGM, "lcv")
end
function updatelin(bhatXl::BallTreeDensity, bXl::BallTreeDensity; N=75)
  dummy = kde!(rand(bhatXl.bt.dims,N),[1.0])
  pGM, = prodAppxMSGibbsS(dummy, [bhatXl, bXl], Union{}, Union{}, 5)
  # error("update -- Not ready")
  return kde!(pGM, "lcv")
end
