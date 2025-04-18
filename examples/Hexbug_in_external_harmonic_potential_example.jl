include(joinpath("..","src","Engine.jl"))


function simulation()

    external_forces = [external_harmonic_force(1,0.5e-1), ABP_3d_propulsion_force(1), self_align_with_v_force(1,1), external_friction_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1])]

    pair_forces = []

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=50
    L = 5.
    initial_state = [ Hexbug([i],[1], [1], [1], [0], [0.1], [0.001], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    
    size = [L,L,2.];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,1e9);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 1e-3, 1e5, Tplot=1e2,fps=120, plot_functions=(plot_points!, plot_directors!, plot_velocity_vectors!)); 
    return sim

end


sim = simulation()
