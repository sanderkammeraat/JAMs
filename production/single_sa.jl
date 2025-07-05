


include(joinpath("..","src","Engine.jl"))

function simulation()

    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,0.2),ABP_perpendicular_angular_noise(1,[0,0,1]),external_harmonic_force(1,0.4))

    pair_forces = []

    
    
    local_dofevolvers =  [overdamped_pq_evolver(1),overdamped_xvf_evolver(1)]
    global_dofevolvers = []
    

    #Initialize state
    L=3.

    initial_state = [PolarParticle3d([1],[1], [1], [1], [1.], [0.3], [0.01], [0,0,0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0])];

    size = [L,L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []
    field_dofevolvers = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers, false,L);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.01,1000,Tplot=10, seed=2, Tsave = 10, save_folder_path=joinpath("/Volumes/T7_Shield/sa/single", "more_noise","not_all_1","simdata"), save_functions = [save_2d_polar_p!], fps=120, plot_functions=[plot_points!, plot_directors!, plot_velocity_vectors!], plotdim=2);
    return sim
end

sim=simulation()



