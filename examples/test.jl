include("../src/Engine.jl")
using Random, Distributions

function simulation()

    Emag = 5e-1
    d = 1e1
    external_forces = [electrode_force(-Emag,d,1),external_friction_force()]

    k=1e0
    pair_forces = [coulomb_force(k)]

    
    
    dofevolvers =  [inertial_evolver!]

    #Initialize state
    N=200
    L=80.



    R = 1

    
    initial_state = [ ChargedParticle3d(i,1,0.1,R,[1.], [rand(Uniform(-L/4, L/4)),rand(Uniform(-L/4, L/4)),rand(Uniform(-d/4, d/4))],[0.,0.,0.],[0.,0.,0.],[0,0,0]) for i=1:N];
    
    size = [L,L,3*d];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,1e9);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    states = Euler_integrator(system, 1e-3, 100000, 100000, 1e2, 120,[new_plot_sized_points!,plot_velocity_vectors!],false);
    return states

end


states = simulation();


    