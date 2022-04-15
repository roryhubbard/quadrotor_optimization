module QuadrotorTrajectory

using JuMP
using Ipopt
using Plots

export test

gr()

function test()
  model = Model(Ipopt.Optimizer)

  # constants
  n = 20 # number of timesteps
  g = 9.8 # gravitational constant
  m = .027 # mass (kg)
  l = .092 # length (m)
  I = 1/12 * m * l^2

  u_min = -5
  u_max = 5

  @variables(model, begin
    X[1:6, 1:n+1]
    u_min <= U[1:2, 1:n] <= u_max # u1 = thrust, u2 = torque
    T >= 0
  end)

  Δt = T / n # discretization length

  # linearized discrete affine state space model (planar quadrotor)
  A = [1 Δt 0  0     0  0
       0  1 0  0 -g*Δt  0
       0  0 1 Δt     0  0
       0  0 0  0     0  0
       0  0 0  0     1 Δt
       0  0 0  0     0  0]
  B = [   0    0
          0    0
          0    0
       Δt/m    0
          0    0
          0 Δt/I]
  affine_term = [0 0 0 -g*Δt 0 0]'

  xi = [-10 0 0 0 0 0]'
  xf = [  0 0 0 0 0 0]'

  @constraints(model, begin
    dynamics[i=1:n], X[:, i+1] .== A * X[:, i] + B * U[:, i] + affine_term
    tilt[i=1:n+1], -π/2 <= X[5, i] <= π/2
    initial_state, X[:, 1]   .== xi
    final_state,   X[:, end] .== xf
  end)

  # minimize time
  @objective(model, Min, T)

  # println(model)
  optimize!(model)

  println(solution_summary(model))
  println("optimal time: ", value(T), " seconds")

  plot(value.(X[1,:]), label="y")
  plot!(value.(X[5,:]), label="θ")
  plot!(xlabel="time (s)")
end

end
