include(joinpath("..","..","src","Engine.jl"))

function simulation()

    
    
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    pair_forces = (soft_disk_force(1,1.),pairAN_force(1,true,1.3, 0., 0., 0.1), pair_nematic_alignment_force(1,2.5,0.0))

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    L =  10

    nx1 =1. 
    ny1 = 0
    nx2 = 1.
    ny2= 0
    R1=1.
    R2 = 1.

    initial_state = PolarParticle3d[ PolarParticle3d([1],[1], [1], [1], [R1], [0.1], [0.00], [-R1 , 0.,0.],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],[nx1, ny1, 0],[0,0,0],[0,0,0])];
    push!(initial_state, PolarParticle3d([2],[1], [1], [1], [R2], [0.1], [0.00], [R2, 0.,0.],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],[nx2,ny2,0],[0,0,0],[0,0,0]))

    display(L)
    sizes = [L,L,4];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    external_forces = (ABP_perpendicular_angular_noise(1,[0,0,1]),)

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,2.5);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 1e4, Tplot=10,fps=Inf,plot_functions=(plot_disks_nematic_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end


sim = simulation()  

@profview sim = simulation() 

@profview_allocs sim = simulation()  

 