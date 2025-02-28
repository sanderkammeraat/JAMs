include("../src/Engine.jl")

function simulation()

    

    pair_forces = [soft_disk_force(1,1)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    L =  70/2
    #initial_state = PolarParticle3d[ PolarParticle3d(i,1, 1, 1, Rs[i], 0.1, 0.001, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


    X = vcat(transpose([xi for xi in -L/2:2:L/2, yi in -L/2:2:L/2])...)
    Y = vcat([xi for xi in -L/2:2:L/2, yi in -L/2:2:L/2]...)
    N = length(X)
    poly=1e-4
    Rs = rand(Uniform(1-poly, 1+poly),N)

    display(X)
    display(Y)

    initial_state = PolarParticle3d[ PolarParticle3d(i,1, 1, 1, Rs[i], 0.1, 0.001, [X[i] , Y[i],0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i in eachindex(X) ];

    pins  = zeros(length(initial_state),3)
    for i=eachindex(initial_state)
        pins[i,:].= initial_state[i].x
    end

    
    sizes = [L+20,L+20,4];
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]),external_harmonic_pinning_force(1,0.1,0,pins))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, false,2.5);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-1, 500, 5e0, 5e0, 120,(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), 2); 
    return sim

end


sim = simulation()  

make_movie(sim, "/Users/kammeraat/test_JAMS/movies/","square_lattice_v4.mp4",(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!),60,2)