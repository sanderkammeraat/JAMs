include(joinpath("..","..","src","Engine.jl"))


function relaxation()

    #pair_forces = (soft_disk_force([1,2],[1. 1. ; 1. 1.]),)
    pair_forces = (soft_disk_force([1,2],[1. 1. ; 1. 1.]),)
    

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver([1,2]),overdamped_pq_xyc_evolver([1,2]))
    global_dofevolvers = []
    field_dofevolvers = []

    Nconf=3000
    Nact = 10
    N = Nconf+Nact
    ϕ = 1.0
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)

    types = vcat(ones(Int64,Nconf), ones(Int64,Nact) .*2)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[types[i]], [1], [1], [Rs[i]], [0.5], [0.001], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    external_forces = []
    display(L)
    sizes = [L,L,4];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []
    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 60, Tplot=nothing,fps=60,plot_functions=(plot_disks_type!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end
rx = relaxation()
function simulation(relaxation)

    

    #pair_forces = (soft_disk_force([1,2],[0.01 1. ; 1. 1.]),)
    pair_forces = (soft_disk_force([1,2],[0.01 1. ; 1. 1.]),)  
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver([1,2]),overdamped_pq_xyc_evolver([1,2]))
    global_dofevolvers = []
    field_dofevolvers = []

    initial_state = relaxation.final_particle_state
    sizes = relaxation.system.sizes
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    external_forces = (ABP_3d_propulsion_force(2),ABP_perpendicular_angular_noise(2,[0,0,1]))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,relaxation.system.rcut_pair_global);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 500, Tplot=10,fps=60,plot_functions=(plot_disks_type!,),plotdim=2)# plot_velocity_vectors!), plotdim=2 )#, record_folder_path=joinpath(homedir(),"soft_gel_05-01-2026"), res= (1000,1000)); 
    return sim;

end


sim = simulation(rx)  


 