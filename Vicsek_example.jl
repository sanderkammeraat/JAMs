include("Engine.jl")
include("LivePlottingFunctions.jl")
using Random, Distributions
function simulation()

    #Initialize state
    N=300
    L=10.
    external_forces = (ABP_2d_propulsion_force(), ABP_2d_angular_noise())

    pair_forces = [ Vicsek_align_force(1)]

    
    
    dofevolvers = [overdamped_evolver!]

    initial_state = [ VicsekParticle(i,0.03,0.1,[rand(Uniform(0, L)) ,rand(Uniform(0,L))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],[0],[0.],1.0,[0.,0.],[0.,0.]) for i=1:N];

    size = [L,L]

    system = System(size, initial_state, external_forces, pair_forces,  dofevolvers, true);


    #%%
    states = Euler_integrator(system, 1, 10000, 1000000, 10, plot_directors!);
    0
end

simulation()

