include(joinpath(pwd(),"src","Engine.jl"))

function simulation()



    #dofevolvers = [inertial_evolver!]
    local_dofevolvers =(overdamped_xvf_evolver([1]),overdamped_pq_evolver([1]))
    global_dofevolvers = []
    L=20
    v0 =  0.0
    Dr=0.001
    initial_state = PolarParticle3d[PolarParticle3d([1],[1], [1], [1], [1],[v0], [Dr], [0 , 0, 0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize(rand(Normal(0, 1),3)).*[1,1,0],[0,0,0],[0,0,0])];

    push!(initial_state, PolarParticle3d([2],[1], [1], [1], [1],[v0], [Dr], [0 , 2, 0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize(rand(Normal(0, 1),3)).*[1,1,0],[0,0,0],[0,0,0]))

    push!(initial_state, PolarParticle3d([3],[1], [1], [1], [1],[v0], [Dr], [2 , 2, 0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize(rand(Normal(0, 1),3)).*[1,1,0],[0,0,0],[0,0,0]))

    push!(initial_state, PolarParticle3d([4],[1], [1], [1], [1],[v0], [Dr], [3 , 1, 0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize(rand(Normal(0, 1),3)).*[1,1,0],[0,0,0],[0,0,0]))

    push!(initial_state, PolarParticle3d([5],[1], [1], [1], [1],[v0], [Dr], [2 , 0,0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize(rand(Normal(0, 1),3)).*[1,1,0],[0,0,0],[0,0,0]))
    


    k_network = [ 0 1 0 0 1 ;
                  0 0 1 0 0; 
                  0 0 0 1 1;
                  0 0 0 0 1;
                  0 0 0 0 0 ]*1

    k_network = k_network + transpose(k_network)

    l_network = [ 0 2 0 0       2;
                  0 0 2 0       0; 
                  0 0 0 sqrt(2) 2;
                  0 0 0 0       sqrt(2);
                  0 0 0 0       0]*1

    l_network = l_network + transpose(l_network)
    sizes = [L,L,4];
    initial_field_state=[]
    field_forces = []
    field_dofevolvers = []
    field_updaters = []

    external_forces = (ABP_3d_propulsion_force(1),ABP_perpendicular_angular_noise(1,[0,0,1]), self_align_with_v_unit_force(1,0.3))

    pair_forces =(spring_network_2d_force(1,l_network, k_network),)


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,10.);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-2, 1e5, Tplot=5e1, fps=120, plot_functions=(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


sim = simulation();

