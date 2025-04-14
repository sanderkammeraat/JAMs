include(joinpath("..","src","Engine.jl"))
using Random, Distributions

function simulation()

    dofevolvers = [overdamped_evolver!]

    #Initialize state
    N=500
    L=4.
    external_forces =  (ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1))


    pair_forces =  (swarm_pos_force(1,1/N, 1),swarm_angular_force(1,1/N,-0.1))
    
    #pair_forces = [swarm_pos_force(1/N, 1)]
    initial_state = [ Swarmalator([i],[1],[0.0],[0.00],rand(Uniform(-L/3, L/3),2),[0.,0.],[0.,0.],[0.,0.],[1.],[1.],[rand(Uniform(-pi, pi))],[0.],[rand(Uniform(-pi, pi))],[0.],[1.0],[0.,0.],[0.,0.],[0,0]) for i=1:N];

    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,L);

    sim = Euler_integrator(system, 0.1, 1e5, Tplot= 5,fps=120, plot_functions=[plot_Swarmalators!]);
    return sim
end

sim=simulation()