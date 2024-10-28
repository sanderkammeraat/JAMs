include("../src/Engine.jl")
using Random, Distributions

function simulation()

    dofevolvers = [overdamped_evolver!]

    #Initialize state
    N=500
    L=4.
    external_forces =  [ABP_2d_propulsion_force()]


    pair_forces =  (swarm_pos_force(1/N, 1),swarm_angular_force(1/N,-0.1))
    
    #pair_forces = [swarm_pos_force(1/N, 1)]
    initial_state = [ Swarmalator(i,0.00,0,[rand(Uniform(L/3, 2*L/3)) ,rand(Uniform(L/3, 2*L/3))],[0.,0.],[0.,0.],1,1,[rand(Uniform(-pi, pi))],[0.],1.0,[0.,0.],[0.,0.]) for i=1:N];

    size = [L,L]

    system = System(size, initial_state, external_forces, pair_forces,  dofevolvers, false);

    states = Euler_integrator(system, 0.1, 1000, 1000000, 10, plot_Swarmalators!);
    0
end

simulation()

