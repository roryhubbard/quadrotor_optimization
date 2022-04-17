module QuadrotorTrajectory

using JuMP, AmplNLWriter, Bonmin_jll, Ipopt, Gurobi
using Plots
gr()

include("../src/BranchAndBound.jl")
using .BranchAndBound

export minimum_time_trajectory

# constants
g = 9.8 # gravity
# crazyflie 2.1 specs
m = .027 # mass (kg)
l = .092 # length (m)
I = 1/12 * m * l^2

struct PlanarQuadrotorLinearizedSS
  A::Matrix
  B::Matrix
  affine_term::Matrix
  function PlanarQuadrotorLinearizedSS(Δt, m::Float64, I::Float64)
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
end

function step!(s::PlanarQuadrotor, u::Vector{Float64}, Δt::Float64)
end

function minimum_time_trajectory()
  #model = Model(Ipopt.Optimizer)
  model = Model(Gurobi.Optimizer)
  #model = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe))

  n = 40 # number of timesteps

  max_control = 10
  tilt_sets = [-π/2 0
                0   π/4
                π/4 π/2]
  M = π/2

  @variables(model, begin
    X[1:6, 1:n+1] # y, ydot, z, zdot, θ, θdot
    -max_control <= U[1:2, 1:n] <= max_control # thrust, torque
    #T >= 0
    z[1:3, 1:n+1], Bin
  end)

  #Δt = T / n # discretization length
  Δt = 0.1

  ss = PlanarQuadrotorLinearizedSS(Δt, m, I)
  A = ss.A
  B = ss.B
  affine_term = ss.affine_term

  xi = [10 0 0 0 0 0]'
  xf = [0  0 0 0 0 0]'

  @constraints(model, begin
    dynamics[i=1:n], X[:, i+1] .== A * X[:, i] + B * U[:, i] + affine_term
    initial_state,     X[:, 1] .== xi
    final_state,     X[:, end] .== xf
    max_visits, sum(z[3, :]) <= 4
  end)

  for i in 1:n+1
    @constraints(model, begin
      X[5, i] >= tilt_sets[1, 1] - M * (1 - z[1, i])
      X[5, i] <= tilt_sets[1, 2] + M * (1 - z[1, i])
      X[5, i] >= tilt_sets[2, 1] - M * (1 - z[2, i])
      X[5, i] <= tilt_sets[2, 2] + M * (1 - z[2, i])
      X[5, i] >= tilt_sets[3, 1] - M * (1 - z[3, i])
      X[5, i] <= tilt_sets[3, 2] + M * (1 - z[3, i])
      sum(z[:, i]) == 1
    end)
  end

  # minimize time
  #@objective(model, Min, T)
  @objective(model, Min, X[1, :]' * X[1, :])

  #println(model)
  optimize!(model)

  println(solution_summary(model))
  #println("optimal time: ", value(T), " seconds")

  plot( value.(X[1,:]), label="y")
  plot!(value.(X[5,:]), label="θ")
  #plot!(value.(U[1,:]), label="thrust")
  #plot!(value.(U[2,:]), label="torque")
  #plot!(xlabel="time (s)")
end

end
