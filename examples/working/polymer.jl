#using Distributed

# n=8

# addprocs(n)

#@everything
include(joinpath("..","..","src","Engine.jl"))

include(joinpath("..","..","io","InitialPositionGenerators.jl"))

function simulation(p,N)


    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    kpar = -1.
    kper= 0.
    f_eq_stretch_force = .7
    pair_forces = (polymer_exterior_soft_disk_force(1,1.),polymer_harmonic_bend_force(1, .3), polymer_harmonic_stretch_force(1,1.,f_eq_stretch_force),polymer_align_director_tangent_force(1,10), polymer_pairAN_force(1,true, true,1.5, kpar, kper, p))
    external_forces =(thermal_translational_noise(1, [0.001, 0.001,0]),)#, ABP_3d_propulsion_force(1))

    

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    pf = 1.
    R = 1
    N_in_pol = N

    Npols = 100
    x, y, radii, pol_ids, ids_in_pol, L = stacked_polymers_at_angle(N_in_pol, Npols, R, pf, f_eq_stretch_force)

    #PolarPolymerParticle3d id type pol_id id_in_pol pol_N
    initial_state = PolarPolymerParticle3d[PolarPolymerParticle3d([id],[1],[pol_ids[id]],[ids_in_pol[id]],[N_in_pol], [1], [1], [radii[id]], [0.3], [0.01], [x[id] , y[id],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for id=1:Npols*N_in_pol];

    display(Npols*N_in_pol)
    display(L)

    sizes = [L,L,2];
    display(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,6.);

    save_folder = "/run/media/martin/HENKESGRFAT/martin/sim_data/p_$p/"
    sim = Euler_integrator(system,0.025, 5000, fps=30, Tplot=nothing, plot_functions=(plot_polymers!, plot_nematic_directors!, plot_velocity_vectors!), plotdim=2, Tsave=40, save_functions=(save_2d_polymer_polar_p!,),save_folder_path = save_folder); 
    return sim;

end 


#@sync @distributed

for p in [.3, 0.4]
    display(p)
    display(Threads.nthreads())
    sim = simulation(p, 10)
    display("Done")
    
end




#sim = simulation(.1, 10) 

# @profview sim = simulation() 

# @profview_allocs sim = simulation()  


 