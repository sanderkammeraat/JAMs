include(joinpath("..","..","src","Engine.jl"))

function simulation()


    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,true,1.3, -1, 0, 0.3), pair_nematic_alignment_force(1,2.5,0.3))


    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = ()
    field_dofevolvers = ()

    N=5000
    ϕ = 1.0
    poly=15e-6
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.1], [0.01], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


    display(L)
    sizes = [L,L,4];
    print(sizes)
    initial_field_state= ()
    field_forces = ()
    field_updaters = ()

    #β=-1 interesting!
    external_forces = () #(ABP_perpendicular_angular_noise(1,[0,0,1]),)

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,3.);

    #Run integrationov
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plot}ting
    sim = Euler_integrator(system,0.05, 1e4,fps=30,Tplot=10,plot_functions=(plot_transparant_disks!,plot_nematic_directors! ),plotdim=2)#, plot_nematic_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end

sim = simulation()  

@profview_allocs simulation()

