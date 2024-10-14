include("Engine.jl")
#include("Forces.jl")
include("Particles.jl")
include("DOFevolvers.jl")
include("LivePlottingFunctions.jl")

using Random, Distributions
using JET
#using BenchmarkTools
function simulation()

    external_forces = []

    #push!(forces, Force("ABP propulsion", "external", Dict("nothing"=>0.),contribute_2d_ABP_propulsion_force!))
    push!(external_forces, ABP_2d_propulsion_force())
    push!(external_forces, ABP_2d_angular_noise())

    pair_forces = []
    

    push!(pair_forces, soft_disk_force())
    print(pair_forces)
    
    dofevolvers = [overdamped_x_evolver!,overdamped_v_evolver!, overdamped_f_evolver!, overdamped_θ_evolver!, overdamped_ω_evolver!]

    #Initialize state
    N=1000
    L=80.
    initial_state = [ PolarParticle2d(i,1,0.3,0.01,[rand(Uniform(0, L)) ,rand(Uniform(0, L))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1.0,1.,1.,[0.,0.],[0.,0.]) for i=1:N];


    size = [L,L];

    system = System(size, initial_state, external_forces, pair_forces , dofevolvers, true);

    #Run integration
    states = Euler_integrator(system, 0.1, 1000, 1000, 10, plot_points!);
    0

end

simulation()

@time simulation()

@profview simulation()

@report_opt simulation()

@benchmark simulation()