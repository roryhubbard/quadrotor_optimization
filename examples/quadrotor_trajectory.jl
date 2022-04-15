module QuadrotorTrajectory

using JuMP
using Ipopt
using Plots
gr()

export minimum_time_trajectory

function minimum_time_trajectory()
  model = Model(Ipopt.Optimizer)

  n = 10 # number of timesteps
  g = 9.8 # gravitational constant

  # crazyflie 2.1 specs
  m = .027 # mass (kg)
  l = .092 # length (m)
  I = 1/12 * m * l^2

  u_min = -5
  u_max = 5
  max_tilt = π/2

  @variables(model, begin
    X[1:6, 1:n+1] # y, ydot, z, zdot, θ, θdot
    u_min <= U[1:2, 1:n] <= u_max # thrust, torque
    T >= 0
  end)

  Δt = T / n # discretization length

  # linearized discrete affine state space model (planar quadrotor)
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

  xi = [-10 0 0 0 0 0]'
  xf = [ 0  0 0 0 0 0]'

  @constraints(model, begin
    dynamics[i=1:n], X[:, i+1] .== A * X[:, i] + B * U[:, i] + affine_term
    tilt[i=1:n+1],   -max_tilt  <= X[5, i] <= max_tilt
    initial_state,     X[:, 1] .== xi
    final_state,     X[:, end] .== xf
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
