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
  X::Vector{Float64} = [y, ydot, z, zdot, θ, θdot]
end

function step!(pq::PlanarQuadrotor, u::Vector{Float64}, Δt::Float64)
  ss = PlanarQuadrotorLinearizedSS(Δt::Float64, m::Float64=m, I::Float64=I)
  pq.X = ss.A * pq.X + ss.B * u + ss.affine_term
  pq.y = pq.X[1]
  pq.ydot = pq.X[2]
  pq.z = pq.X[3]
  pq.zdot = pq.X[4]
  pq.θ = pq.X[5]
  pq.θdot = pq.X[6]
end

function minimum_time_trajectory()
  #model = Model(Ipopt.Optimizer)
  model = Model(Gurobi.Optimizer)
  #model = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe))

  n = 40 # number of timesteps
  max_control = 10
  tilt_sets = [-π/4 0
                0   π/8
                π/8 π/4]
  M = π/2

  @variables(model, begin
    #X[1:6, 1:n+1] # y, ydot, z, zdot, θ, θdot
    y[1:n+1]
    ydot[1:n+1]
    -1 <= z[1:n+1] <= 1
    -1 <= zdot[1:n+1] <= 1
    θ[1:n+1]
    θdot[1:n+1]
    -max_control <= U[1:2, 1:n] <= max_control # thrust, torque
    #T >= 0
    bv[1:3, 1:n+1], Bin
  end)

  X = vcat(y', ydot', z', zdot', θ', θdot') # y, ydot, z, zdot, θ, θdot

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
    max_visits, sum(bv[3, :]) <= 4
  end)

  for i in 1:n+1
    @constraints(model, begin
      X[5, i] >= tilt_sets[1, 1] - M * (1 - bv[1, i])
      X[5, i] <= tilt_sets[1, 2] + M * (1 - bv[1, i])
      X[5, i] >= tilt_sets[2, 1] - M * (1 - bv[2, i])
      X[5, i] <= tilt_sets[2, 2] + M * (1 - bv[2, i])
      X[5, i] >= tilt_sets[3, 1] - M * (1 - bv[3, i])
      X[5, i] <= tilt_sets[3, 2] + M * (1 - bv[3, i])
      sum(bv[:, i]) == 1
    end)
  end

  # minimize time
  #@objective(model, Min, T)
  @objective(model, Min, sum(X[1, :].^2))

  #println(model)
  optimize!(model)

  println(solution_summary(model))
  #println("optimal time: ", value(T), " seconds")

  p1 = plot( value.(X[1,:]), label="y")
  plot!(p1, value.(X[2,:]), label="ydot")
  plot!(p1, value.(X[3,:]), label="z")
  plot!(p1, value.(X[4,:]), label="zdot")
  plot!(p1, value.(X[5,:]), label="θ")
  plot!(p1, value.(X[6,:]), label="θdot")

  p2 = plot(value.(U[1,:]), label="thrust")
  plot!(p2, value.(U[2,:]), label="torque")

  plot(p1, p2, layout=(2,1))
  #plot!(value.(U[1,:]), label="thrust")
  #plot!(value.(U[2,:]), label="torque")
  #plot!(xlabel="time (s)")

#  anim = @animate for i = 1:num_time_steps
#    plot(q[1, :], q[2, :], xlim=(-1.1, 1.1), ylim=(-1.1, 1.1))
#    plot!([q[1, i]], [q[2, i]], marker=(:hex, 6))
#  end
#
#  gif(anim, "img/mpc1.gif", fps = 30)
end

end
