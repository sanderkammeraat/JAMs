include("../src/Engine.jl")

function simulation()

    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,0.4),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force(1,1), chain_force(1,1e-1,2)]

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=20
    L=4 * N*1.
    v0 =  0.3
    Dr=0.01
    initial_state = [PolarParticle3d(1,1, 1, 1, 1,v0, Dr, [0 , -1,0],[0,0,0], [0,0,0],[0,0,0],normalize(rand(Normal(0, 1),3)).*[1,1,0],[0,0,0],[0,0,0])];

    for i in 2:N
        
        push!(initial_state, PolarParticle3d(i,1, 1, 1, 1,v0, Dr, [0 , -2*i+1,0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
        print(-2*i+1)
    end

    size = [2.5*L,2.5*L,4];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,10.);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-1, 1e5, 1e10, 1e0, 120,(plot_sized_points!, plot_directors!, plot_velocity_vectors!),2); 
    return sim

end


sim = simulation();

