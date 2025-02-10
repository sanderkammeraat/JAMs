include("../src/Engine.jl")
using Random, Distributions

function simulation()

    external_forces = [ABP_2d_angular_noise(1)]

    pair_forces = [soft_disk_force(1,1)]

    field_forces = [field_propulsion_force(1,1e-2,0.01)]
    

    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1000
    Lx=300.
    Ly=300.
    L=min(Lx,Ly)
    Li = 50
    poly = 0.0000002


    initial_particle_state = [ PolarParticle2d(i,1,1,0.0,0.0001,[rand(Uniform(-Li, Li)) ,rand(Uniform(-Li,Li))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.]) for i=1:N];
    
    size = [Lx,Ly];

    lbin = 1
    print(lbin)

    x_bin_centers = range(start=-Lx/2, stop=Lx/2, step=lbin).+lbin/2
    y_bin_centers = range(start=-Ly/2, stop=Ly/2, step=lbin).+lbin/2
    bin_centers = [x_bin_centers, y_bin_centers]
    print(x_bin_centers)

    v00=0.4
    field_updaters = [PeriodicDiffusion([1],2e-2), EdgeSet([1],v00)]

    C = ones(length(x_bin_centers), length(y_bin_centers))*v00
    
    initial_field_state=[GeneralField2d(1,bin_centers,C, C.*0, C.*0)]

    system = System(size, initial_particle_state,initial_field_state,external_forces, pair_forces, field_forces, field_updaters,dofevolvers, true, 1e1);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    plot_functions = (plot_sized_points!, plot_directors!, plot_velocity_vectors!,plot_field_magnitude!)
    sim = Euler_integrator(system, 0.1, 1e7, 1e7,5, 120, plot_functions,2);
    #particle_states,field_states = Euler_integrator(system, 0.1, 10, 1e8,0, 0, plot_functions,false);
    return sim
end


sim = simulation()