include(joinpath("..","src","Engine.jl"))

function simulation(v0, Dr, J, Tplot, Tsave,tend)

    

    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    pair_forces =[]#(soft_disk_force(1,0),)


    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    #Number of particles
    N=1
    ϕ = 1.0
    poly=15e-6
    Rs = rand(Uniform(1-poly, 1+poly),N)

    L =  5

    sigma = 1
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [v0], [Dr], [-1+sigma*rand(Normal(0, 1)) , sigma*rand(Normal(0, 1)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    xa = -1
    xb = 1
    ka = 2
    kb = 2


    initial_field_state=[]

    sizes = [L,L,4];
    #initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    external_forces = (external_double_gaussian_force(1,ka,kb, [xa,0,0], [xb,0,0]),ABP_perpendicular_angular_noise(1,[0,0,1]),self_align_with_v_unit_force(1,J),ABP_3d_propulsion_force(1))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,3.);

    #Run integrationov
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting

    #save_folder = "/Users/kammeraat/dwsa/single/simdata/v0_$v0/Dr_$Dr/J_$J/"
    save_folder = "/Users/kammeraat/sa_double_well/"
    sim = Euler_integrator(system,0.01,tend, Tsave=Tsave, fps=60,Tplot=Tplot,plot_functions=(plot_points!, plot_velocity_vectors!), plotdim=2, save_folder_path = save_folder, save_functions = (save_2d_polar_p!,)); 
    return sim;

end

#v0, Dr, J, Tplot, Tsave

v0 = 0.6
Dr = 0.01
J = 1
Tplot = 1 #plot every timestep, set to Tplot=nothing to turn off plotting

Tsave = 100 #save every nth timestep, set to Tsave=nothing to turn off saving

tend = 1e3 #for how long to run the simulation (in units of simulation time)

simulation(v0, Dr ,J, Tplot, Tsave,tend)
