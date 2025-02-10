include("../src/Engine.jl")

function simulation()

    external_forces = [ ABP_3d_propulsion_force(1), self_align_with_v_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]),external_harmonic_force(1,2e-3)]

    pair_forces = [soft_disk_force(1,1)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=1000
    ϕ = 0.3
    L =  sqrt(N *  π * 1^2 / ϕ)
    initial_state = [ PolarParticle3d(i,1, 1, 1, 1, 0.3, 0.01, [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0]) for i=1:N ];

    
    size = [L,L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,1e1);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-1, 1e5, 1e10, 5e0, 120,(plot_sized_points!, plot_directors!, plot_velocity_vectors!), 2); 
    return sim

end


sim = simulation()

