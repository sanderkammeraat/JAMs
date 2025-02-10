include("../src/Engine.jl")
using Random, Distributions
function simulation()

    #Initialize state
    N=3000
    L=32.
    external_forces = (ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1))

    pair_forces = [ Vicsek_align_force(1,1)]

    
    
    dofevolvers = [overdamped_evolver!]

    initial_state = [ VicsekParticle(i,1,0.5,0.05,[rand(Uniform(-L/2, L/2)) ,rand(Uniform(-L/2,L/2))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],[0],[0.],1.0,[0.,0.],[0.,0.]) for i=1:N];

    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,3e0);


    #%%
    sim = Euler_integrator(system, 1, 10000, 1000000, 1,120, [plot_directors!]);
    return sim
end

sim=simulation()

