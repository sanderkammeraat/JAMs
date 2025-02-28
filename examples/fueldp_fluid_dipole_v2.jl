include("../src/Engine.jl")

function simulation()

    external_forces = (ABP_3d_propulsion_force(1),self_align_with_v_force(1,0.3),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = (soft_disk_force(1,1),fluid_dipole_3d_force(1,1e-1,6))

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=1000
    ϕ = 0.2
    L =  sqrt(N *  π * 1^2 / ϕ)
    initial_state = PolarParticle3d[ PolarParticle3d(i,1, 1, 1, 1, 0.3, 0.001, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    
    size = [L,L,2];
    lbin = 1
    print(lbin)
    Lx=L
    Ly=L
    x_bin_centers = range(start=-Lx/2, stop=Lx/2, step=lbin).+lbin/2
    y_bin_centers = range(start=-Ly/2, stop=Ly/2, step=lbin).+lbin/2
    z_bin_centers = [0.]
    bin_centers = [x_bin_centers, y_bin_centers,z_bin_centers]
    print(x_bin_centers)

    v00=0.3

    C = ones(length(x_bin_centers), length(y_bin_centers))*v00
    
    initial_field_state=[FuelField2d(1,bin_centers,C, C.*0, C.*0)]
    field_forces =[]# [field_propulsion_3d_force(1,1e-2,0.01)]
    field_updaters = [PeriodicDiffusion(1,8e-3)]




    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,20.);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,5e-2, 2e2, 2e1, 2e1, 120,(plot_field_magnitude!,plot_disks!, plot_director_points!, plot_velocity_vectors!), 2); 
    return sim

end


sim = simulation()

make_movie(sim, "/Users/kammeraat/test_JAMS/movies/","fueldp_fluid_dipole_tstop_1000.mp4",(plot_disks!, plot_directors!, plot_velocity_vectors!,plot_field_magnitude!),60,2)