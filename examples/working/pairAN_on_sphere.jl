include(joinpath("..","..","src","Engine.jl"))
include(joinpath("..","..","io","InitialPositionGenerators.jl"))



function simulation()


    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3

    kpar = -1
    kper=0
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,true,1.3, kpar, kper, 0.1), pair_nematic_alignment_force(1,2.5,0.1))


    #dofevolvers = [inertial_evolver!]
    global_dofevolvers = ()
    field_dofevolvers = ()

    N=3000
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    R = 20
    local_dofevolvers = (overdamped_xvf_Rc_evolver(1,R),overdamped_pq_Rc_evolver(1,R))

    L =  3*R

    positions, polarities = random_on_sphere(R,N)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.2], [0.000], positions[i],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],polarities[i],[0,0,0],[0,0,0]) for i=1:N ];


    display(L)
    sizes = [L,L,L];
    print(sizes)
    initial_field_state= ()
    field_forces = ()
    field_updaters = ()

    #β=-1 interesting!
    external_forces = (ABP_3d_propulsion_force(1),ABP_3d_angular_noise(1))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,3.);

    #Run integrationov
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plot}ting
    sim = Euler_integrator(system,0.01, 1e4,fps=60,Tplot=10,plot_functions=(plot_sphere!,plot_sized_points!,plot_nematic_directors! ),plotdim=3)#, plot_nematic_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end

sim = simulation()  

@profview simulation()



