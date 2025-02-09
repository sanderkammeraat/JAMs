include("../src/Engine.jl")
using Random, Distributions

function simulation()

    dofevolvers = [overdamped_evolver!]

    #Initialize state
    N=500
    L=4.
    external_forces =  (ABP_2d_propulsion_force(), ABP_2d_angular_noise())


    pair_forces =  (swarm_pos_force(1/N, 1),swarm_angular_force(1/N,-0.1))
    
    #pair_forces = [swarm_pos_force(1/N, 1)]
    initial_state = [ Swarmalator(i,0.0,0.00,rand(Uniform(-L/3, L/3),2),[0.,0.],[0.,0.],1.,1.,[rand(Uniform(-pi, pi))],[0.],[rand(Uniform(-pi, pi))],[0.],1.0,[0.,0.],[0.,0.]) for i=1:N];

    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,1e9);

    sim = Euler_integrator(system, 0.1, 1000, 1000000, 5,120, [plot_Swarmalators!]);
    return sim
end

sim=simulation()

