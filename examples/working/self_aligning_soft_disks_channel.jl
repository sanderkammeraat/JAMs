include(joinpath("..","..","src","Engine.jl"))

function simulation()

    

    pair_forces = (soft_disk_force(1,0.4),pair_polar_alignment_force(1,3,1.))

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    N=1000
    ϕ = 0.9
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.1], [0.01], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


    display(L)
    sizes = [L,L,4];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,-0.5),ABP_perpendicular_angular_noise(1,[0,0,1]))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,5.);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 1e5, Tplot=10,fps=120,plot_functions=(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end


sim = simulation()  

@profview sim = simulation() 

@profview_allocs sim = simulation()  

 