include(joinpath("..","..","src","Engine.jl"))

function simulation()

    

    pair_forces = (soft_disk_force(1,1.),)
    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,0.1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    N=10000
    ϕ = 0.2
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)
    Ly = L/2
    Lx = Ly*2

    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.3], [0.001], [rand(Uniform(-Lx/2, Lx/2)) , rand(Uniform(-Ly/2,Ly/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


    display(L)
    sizes = [Lx,Ly,1];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 100, Tplot=10,fps=60,plot_functions=(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2, Tsave=nothing, save_folder_path = "/Users/kammeraat/test_saving_speed/", save_functions=[save_2d_polar_p!]); 
    return sim;

end


sim = simulation()  

@profview sim = simulation() 

@profview_allocs sim = simulation()  

 