#using Distributed

# n=8

# addprocs(n)

#@everything
include(joinpath("..","..","src","Engine.jl"))

include(joinpath("..","..","io","InitialPositionGenerators.jl"))


sim_folder_name = "sim_data"

#for windows
path_data = joinpath("E:", "martin", sim_folder_name)

#for linux
#path_data = joinpath("/run/media/martin/HENKESGRFAT/martin", sim_folder_name)

function simulation(p, N_in_pol, kpar, kbend)


    f_eq_stretch_force = .75
    krep = 1.

    pair_forces = (soft_disk_force(1,krep), polymer_harmonic_bend_force(1, kbend), polymer_pairAN_force(1,false, false, false, 1.5, kpar, 0., p), polymer_harmonic_stretch_force(1,krep, f_eq_stretch_force))
    #pair_forces = (soft_disk_force(1,krep), polymer_harmonic_stretch_force(1,krep, f_eq_stretch_force), polymer_harmonic_bend_force(1, kbend), polymer_pair_polar_nematic_force(1, false, 1.5, p))

    external_forces = (thermal_translational_noise(1, [0.001, 0.001, 0]),)#, ABP_3d_propulsion_force(1))


    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),)
    global_dofevolvers = (polymer_p_set_evolver(1),)
    field_dofevolvers = ()

    pf = 1.
    R = 1
    numb_particles = 400

    Npols = convert(Int64, floor(numb_particles/N_in_pol))
    x, y, radii, pol_ids, ids_in_pol, L = stacked_polymers_at_angle(N_in_pol, Npols, R, pf, f_eq_stretch_force, random_polarity = false)

    #PolarPolymerParticle3d id type pol_id id_in_pol pol_N #####,overdamped_pq_xyc_evolver(1)
    #initial_state = PolarPolymerParticle3d[PolarPolymerParticle3d([id],[1],[pol_ids[id]],[ids_in_pol[id]],[N_in_pol], [1], [1], [radii[id]], [0.], [0.01], [x[id] , y[id],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for id=1:Npols*N_in_pol];
    
    file = jldopen(raw"E:\martin\sim_data\p_0.4,N_10\JAMs_final_state.jld2", "r")
    older_simulation = file["SIM"]
    close(file)

    initial_state = older_simulation.final_particle_state
    
    display(Npols*N_in_pol)

    sizes = [L,L,2];
    print(sizes)
    initial_field_state=()
    field_forces = ()
    field_updaters = ()

    
    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,6.);

    save_folder = joinpath(path_data, "p_$p,N_$N_in_pol/")
    sim = Euler_integrator(system, 0.01, 1000, fps=30, Tplot=40, plot_functions=(plot_polymers!, plot_nematic_directors!), plotdim=2, Tsave=nothing, save_functions=(save_2d_polymer_polar_p!,), save_folder_path = save_folder);#, record_folder_path = pwd()); 
    return sim;

end 


#@sync @distributed


#p, N_in_pol, kpar, kbend
sim = simulation(.4,  10, -1., 3.) 

# @profview sim = simulation() 

# @profview_allocs sim = simulation()  


 