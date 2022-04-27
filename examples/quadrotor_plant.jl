module QuadrotorPlant

export DynamicsLinearized

# constants
g = 9.8 # gravity
# crazyflie 2.1 specs
m = .027 # mass (kg)
l = .092 # length (m)
I = 1/12 * m * l^2

struct DynamicsNonlinear
  function DynamicsNonlinear(Δt, m::Float64, I::Float64)
  end
end

struct DynamicsLinearized
  A::Matrix
  B::Matrix
  affine_term::Matrix
  function DynamicsLinearized(Δt, m::Float64, I::Float64)
    A = [1  Δt 0  0   0    0
         0  1  0  0  -g*Δt 0
         0  0  1  Δt  0    0
         0  0  0  1   0    0
         0  0  0  0   1    Δt
         0  0  0  0   0    1]
    B = [0    0
         0    0
         0    0
         Δt/m 0
         0    0
         0    Δt/I]
    affine_term = [0 0 0 -g*Δt 0 0]'
    new(A, B, affine_term)
  end
end

mutable struct PlanarQuadrotor
  y::Float64
  ydot::Float64
  z::Float64
  zdot::Float64
  θ::Float64
  θdot::Float64
  function PlanarQuadrotor(Δt::Float64, m::Float64, I::Float64, linearize=false)
  end
end

function step!(pq::PlanarQuadrotor, u::Vector{Float64}, Δt::Float64)
  ss = DynamicsLinearized(Δt::Float64, m::Float64, I::Float64)
  X = ss.A * pq.X + ss.B * u + ss.affine_term
  pq.y    = X[1]
  pq.ydot = X[2]
  pq.z    = X[3]
  pq.zdot = X[4]
  pq.θ    = X[5]
  pq.θdot = X[6]
end

end
