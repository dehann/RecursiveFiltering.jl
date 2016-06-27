
type SysState
    X::Array{Float64,1}
    P::Array{Float64,2}
end

function c2d!(M1::Array{Float64,2},
              M2::Array{Float64,2},
              F::Array{Float64,2},
              G::Array{Float64,2},
              Qc::Array{Float64,2},
              dt::Float64;
              expmfnc::Function=expm)
    # more accurate explicit matrix exponential approach --
    # computationally efficient approaches can be done here to lighten the load

    # M1 = [[F';G']';zeros(12,21)]
    # M1 = zeros(21,21)
    M1[1:9,1:9] = F
    M1[1:9,10:21] = G
    M1[10:21,1:21] = 0.0

    Md1 = expmfnc(M1*dt)
    Phi = Md1[1:9,1:9]
    Gamma = Md1[1:9,10:21]

    #M2 = [[-F';(G*Qc*G')']';[zeros(9,9);F]']
    # M2 = zeros(18,18)
    M2[1:9,1:9] = -F
    M2[1:9,10:18] = (G*Qc*G')
    M2[10:18,10:18] = F'
    M2[10:18,1:9] = 0.0

    Md2 = expmfnc(M2*dt)
    Qd = Phi*(Md2[1:9,10:18])
    return Phi, Gamma, Qd
end


# cpp: ReturnTypePhiGammaQd c2d(const Eigen::MatrixXd &F, const Eigen::MatrixXd &G, const Eigen::MatrixXd &Qc, const double &dt) {
function c2d(F::Array{Float64,2}, G::Array{Float64,2}, Qc::Array{Float64,2}, dt::Float64;
            expmfnc::Function=expm)
    # more accurate explicit matrix exponential approach --
    # computationally efficient approaches can be done here to lighten the load

    # M1 = [[F';G']';zeros(12,21)] # easy concat
    M1 = zeros(21,21)
    # Eigen has M1.topLeftCorner(9,9) = F;
    M1[1:9,1:9] = F
    M1[1:9,10:21] = G

    Md1 = expmfnc(M1*dt)
    # Eigen has Phi = Md1.topLeftCorner(9,9);
    Phi = Md1[1:9,1:9]
    Gamma = Md1[1:9,10:21]

    #M2 = [[-F';(G*Qc*G')']';[zeros(9,9);F]'] # easy concat
    M2 = zeros(18,18)
    M2[1:9,1:9] = -F
    M2[1:9,10:18] = (G*Qc*G')
    M2[10:18,10:18] = F'

    Md2 = expmfnc(M2*dt)
    Qd = Phi*(Md2[1:9,10:18])
    return Phi, Gamma, Qd
end

function filterpropagate!(st::SysState, wQb::Quaternion, dt::Float64, Qcont::Array{Float64,2})
    F, G = FG(wQb)
    Phi, Gamma, Qd = c2d(F, G, Qcont, dt)
    # ultimately move M1 and M2 declaration all the way up to above the iteration loop to allow memory reuse

    # estimated error state
    st.X = Phi*st.X # + Gamma*noiseinput, which we can only approximate
    # noise increases uncertainty (error covariance)
    st.P = Phi*st.P*Phi' + Qd
    nothing
end

function filterupdate!(st::SysState, r::Array{Float64,1}, H::Array{Float64,2}, R::Array{Float64,2})
    # K = P H' (A.inv) = P (H'*Ainv) = P (A' \ H)'
    At = (H*st.P*(H') + R)'
    K = st.P * ((At \ H)')
    # inverse way
    #K = st.P*H'*((H*st.P*(H') + R) \ eye(size(R,1))) # inverse
    st.X += K*r
    KH = K*H
    P = (eye(size(KH,1)) - KH)*st.P
    st.P = P
    nothing
end
