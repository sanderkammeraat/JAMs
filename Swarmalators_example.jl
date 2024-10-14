include("Engine.jl")
#include("Forces.jl")
include("Particles.jl")
include("DOFevolvers.jl")
include("LivePlottingFunctions.jl")
using Random, Distributions
using JET
function simulation()

    dofevolvers = [overdamped_x_evolver!,overdamped_v_evolver!, overdamped_f_evolver!, overdamped_θ_evolver!, overdamped_ω_evolver!]


    #Initialize state
    N=1000
    L=5.
    external_forces = []

    push!(external_forces, ABP_2d_propulsion_force())
    #push!(forces, Force("ABP angular noise", "external", Dict(), contribute_2d_ABP_angular_noise!))

    pair_forces = []
    push!(pair_forces, swarm_pos_force(1/N, 1.) )
    push!(pair_forces, swarm_angular_force(1/N,-0.1) )

    
    initial_state = [ Swarmalator(i,0.01,0,[rand(Uniform(L/3, 2*L/3)) ,rand(Uniform(L/3, 2*L/3))],[0.,0.],[0.,0.],1,1,[rand(Uniform(-pi, pi))],[0.],1.0,[0.,0.],[0.,0.]) for i=1:N];

    size = [L,L]

    system = System(size, initial_state, external_forces, pair_forces,  dofevolvers, false);


    #%%
    states = Euler_integrator(system, 0.1, 1000, 1000, 10, plot_Swarmalators!);
    0
end

simulation()

@time simulation()

@profview simulation()

@report_opt simulation()