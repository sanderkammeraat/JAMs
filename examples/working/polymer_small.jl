include(joinpath("..","..","src","Engine.jl"))

include(joinpath("..","..","io","InitialPositionGenerators.jl"))

function simulation()


    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    kpar = 0
    kper= 1.
    pair_forces = (soft_disk_force(1,1),polymer_harmonic_bend_force(1,0.1), polymer_harmonic_stretch_force(1,1.),polymer_align_director_tangent_force(1,0.001), polymer_pairAN_force(1,true, false,1.3, kpar, kper, 0.3))
    external_forces =[]

    

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = ()
    field_dofevolvers = ()


    L=40

    N_in_pol=10

    x = collect(0:2:2*N_in_pol) .- N_in_pol

    y = zeros(N_in_pol)

    #PolarPolymerParticle3d id type pol_id id_in_pol pol_N
    initial_state = PolarPolymerParticle3d[PolarPolymerParticle3d([id],[1],[1],[id],[N_in_pol], [1], [1], [1], [0.3], [0.01], [x[id] , y[id],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for id=1:N_in_pol];

    display(L)

    sizes = [L,L,2];
    print(sizes)
    initial_field_state=()
    field_forces = ()
    field_updaters = ()

    #β=-1 interesting!
    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,6.);

    sim = Euler_integrator(system,0.01, 1e4, Tplot=10,plot_functions=(plot_polymers!, plot_nematic_directors!, plot_velocity_vectors!), plotdim=2, Tsave=nothing, save_functions=(save_2d_polymer_polar_p!,),save_folder_path = joinpath(homedir(),"ANP","demos","N_in_pol_$(N_in_pol)","k_par_$(kpar)","k_per_$(kper)"),res=(1000,1000)); 
    return sim;

end 

sim = simulation()  

@profview sim = simulation() 

@profview_allocs sim = simulation()  

 