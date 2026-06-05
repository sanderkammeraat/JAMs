include(joinpath("..","..","src","Engine.jl"))
include(joinpath("..","..","io","InitialPositionGenerators.jl"))

function simulation()

    pair_forces = (soft_disk_force(1,1.),)
    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,10),ABP_3d_angular_noise(1))

    #dofevolvers = [inertial_evolver!]
    
    global_dofevolvers = []
    field_dofevolvers = []

    N=500
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    R = 10
    local_dofevolvers = (overdamped_xvf_Rc_evolver(1,R),overdamped_pq_Rc_evolver(1,R))

    L =  2.5*R
    Ly = L
    Lx = L
    Lz = L
    positions, polarities = random_on_sphere(R,N)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.2], [0.000], positions[i],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],polarities[i],[0,0,0],[0,0,0]) for i=1:N ];


    display(L)
    sizes = [Lx,Ly,Lz];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,3*(1+poly));

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.005, 1e4, sbs=false, Tplot=100,fps=60,plot_functions=(plot_sphere!,plot_points!, plot_directors!, plot_velocity_vectors!,plot_trajectories!), Tsave=nothing, save_folder_path = "/Users/kammeraat/test_saving_speed/", save_functions=[save_2d_polar_p!]); 
    return sim;

end
simulation()
