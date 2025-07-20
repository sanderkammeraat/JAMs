include(joinpath(pwd(),"src","Engine.jl"))

function simulation()

    external_forces = (ABP_3d_propulsion_force(1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force(1,1), chain_force(1,1,1), tangential_propulsion_force(1)]

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_evolver(1))
    global_dofevolvers = []
    N=100
    L=4 * N*1.
    v0 =  0.1
    Dr=0.00001
    initial_state = [PolarParticle3d([1],[1], [1], [1], [1],[v0], [Dr], [0 , N ,0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize(rand(Normal(0, 1),3)).*[1,1,0],[0,0,0],[0,0,0])];
    multiplier=1
    for i in 2:N
    
        #multiplier*=-1
        
        push!(initial_state, PolarParticle3d([i],[1], [1], [1], [1],[multiplier*v0], [Dr], [0 ,  N  -2*(i-1),0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
        print(-2*i+1)
    end

    sizes = [L,L,2];
    initial_field_state=[]
    field_forces = []
    field_dofevolvers = []
    field_updaters = []


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,4.);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-2, 1e5, Tplot=5e1, fps=120, plot_functions=(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


sim = simulation();

