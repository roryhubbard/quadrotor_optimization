module QuadrotorTrajectory

using JuMP, AmplNLWriter, Bonmin_jll, Ipopt, Gurobi
using Plots
gr()

include("quadrotor_plant.jl")
using .QuadrotorPlant
#include("../src/BranchAndBound.jl")
#using .BranchAndBound

#model = Model(Ipopt.Optimizer)
#model = Model(() -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe))

# constants
g = 10.0 # gravity
m = 1.0 # mass (kg)
l = 1.0 # length (m)
I = 1.0 # rotational intertia

function maximize_alititude(n::Int, n_max_thrust_limit::Int)
  model = Model(Gurobi.Optimizer)

  Δt = 0.1
  thrust_sets = [ 0  20
                 20  30
                 30  40]
  Mu = 40

  @variables(model, begin
    X[1:6, 1:n+1]
    0 <= U[1:2, 1:n] <= 50
    bv_u[1:3,   1:n], Bin
  end)

  ss = DynamicsLinearized(Δt, m, I)
  A = ss.A
  B = ss.B
  affine_term = ss.affine_term

  xi = [0 0 0 0 0 0]'

  @constraints(model, begin
    dynamics[i=1:n], X[:, i+1] .== A * X[:, i] + B * U[:, i] + affine_term
    initial_state,     xi .<= X[:, 1] .<= xi
    max_visits, sum(bv_u[3, :]) <= n_max_thrust_limit
  end)

  for i in 1:n
    @constraints(model, begin
      U[1, i] >= thrust_sets[1, 1] - Mu * (1 - bv_u[1, i])
      U[1, i] <= thrust_sets[1, 2] + Mu * (1 - bv_u[1, i])
      U[1, i] >= thrust_sets[2, 1] - Mu * (1 - bv_u[2, i])
      U[1, i] <= thrust_sets[2, 2] + Mu * (1 - bv_u[2, i])
      U[1, i] >= thrust_sets[3, 1] - Mu * (1 - bv_u[3, i])
      U[1, i] <= thrust_sets[3, 2] + Mu * (1 - bv_u[3, i])
      sum(bv_u[:, i]) == 1
    end)
  end

  @objective(model, Min, -X[3, end])
  optimize!(model)

  value.(X[3, :]), value.(U[1, :])
end

function maximize_alititude_integer_study()
  n = 10 # number of timesteps
  i = [1, 5, 9]
  a = Vector{Vector{Float64}}()
  b = Vector{Vector{Float64}}()
  for i_ in i
    rx, ru = maximize_alititude(n, i_)
    push!(a, rx)
    push!(b, ru)
  end
  p1 = plot(a[1], label="z₃=1", ylabel="altitude (m)", legend_position=:topleft)
  plot!(p1, a[2], label="z₃=5", ylabel="altitude (m)", legend_position=:topleft)
  plot!(p1, a[3], label="z₃=9", ylabel="altitude (m)", legend_position=:topleft)
  p2 = plot(b, ylabel="thrust (N)", legend_position=false)
  fig = plot(p1, p2, layout=(2, 1), xlabel = "timesteps")
  savefig(fig, "fig")
end

function minimum_time_trajectory()
  model = Model(Gurobi.Optimizer)

  n = 10 # number of timesteps
  Δt = 0.1
  tilt_sets = [-π/4 -π/8
               -π/8  π/8
                π/8  π/4]
  Mθ = π/2

  @variables(model, begin
    X[1:6, 1:n+1]
    0 <= U[1:2, 1:n] <= 50
    bv_θ[1:3, 1:n+1], Bin
  end)

  ss = DynamicsLinearized(Δt, m, I)
  A = ss.A
  B = ss.B
  affine_term = ss.affine_term

  xi = [0 0 0 0 0 0]'

  @constraints(model, begin
    dynamics[i=1:n], X[:, i+1] .== A * X[:, i] + B * U[:, i] + affine_term
    initial_state,     xi .<= X[:, 1] .<= xi
    max_visits, sum(bv_u[3, :]) <= 5
  end)

  for i in 1:n+1
    @constraints(model, begin
      X[5, i] >= tilt_sets[1, 1] - Mθ * (1 - bv[1, i])
      X[5, i] <= tilt_sets[1, 2] + Mθ * (1 - bv[1, i])
      X[5, i] >= tilt_sets[2, 1] - Mθ * (1 - bv[2, i])
      X[5, i] <= tilt_sets[2, 2] + Mθ * (1 - bv[2, i])
      X[5, i] >= tilt_sets[3, 1] - Mθ * (1 - bv[3, i])
      X[5, i] <= tilt_sets[3, 2] + Mθ * (1 - bv[3, i])
      sum(bv[:, i]) == 1
    end)
  end

  # minimize time
  @objective(model, Min, sum(X[3, :].^2))

  optimize!(model)

#  p1 = plot(value.(X[1,:]), label="y")
#  plot!(p1, value.(X[2,:]), label="z")
#  plot!(p1, value.(X[3,:]), label="θ")
#  plot!(p1, value.(X[4,:]), label="ydot")
#  plot!(p1, value.(X[5,:]), label="zdot")
#  plot!(p1, value.(X[6,:]), label="θdot")

#  p2 = plot(value.(U[1,:]), label="thrust")
#  plot!(p2, value.(U[2,:]), label="torque")

#  plot(p1, p2, layout=(2,1))

#  anim = @animate for i = 1:num_time_steps
#    plot(q[1, :], q[2, :], xlim=(-1.1, 1.1), ylim=(-1.1, 1.1))
#    plot!([q[1, i]], [q[2, i]], marker=(:hex, 6))
#  end
#
#  gif(anim, "img/mpc1.gif", fps = 30)
end

end
