#using Distributed

# n=8

# addprocs(n)

#@everything
include(joinpath("..","..","src","Engine.jl"))

include(joinpath("..","..","io","InitialPositionGenerators.jl"))

function simulation(p, kpar, kbend)


    f_eq_stretch_force = .75
    krep = 1.
    pair_forces = (soft_disk_force(1,krep), polymer_harmonic_bend_force(1, kbend), polymer_harmonic_stretch_force(1,krep, f_eq_stretch_force), polymer_pairAN_force(1,true, true, false, 1.5, kpar, 0., p))
    external_forces = (thermal_translational_noise(2, [0.001, 0.001, 0]),)#, ABP_3d_propulsion_force(1))


    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),)
    global_dofevolvers = (polymer_p_set_evolver(1),)
    field_dofevolvers = ()

    pf = 1.
    R = 1
    N_in_pol = 10

    Npols = 151
    x, y, radii, pol_ids, ids_in_pol, L = stacked_polymers_at_angle(N_in_pol, Npols, R, pf, f_eq_stretch_force)

    #PolarPolymerParticle3d id type pol_id id_in_pol pol_N #####,overdamped_pq_xyc_evolver(1)
    initial_state = PolarPolymerParticle3d[PolarPolymerParticle3d([id],[1],[pol_ids[id]],[ids_in_pol[id]],[N_in_pol], [1], [1], [radii[id]], [0.3], [0.01], [x[id] , y[id],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for id=1:Npols*N_in_pol];

    display(Npols*N_in_pol)

    sizes = [L,L,2];
    print(sizes)
    initial_field_state=()
    field_forces = ()
    field_updaters = ()

    
    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,6.);

    save_folder = joinpath(raw"C:\Users\gabri\Documents\Travail-Etude\Master's Theoretical Physics Leiden\Research Project\Data", "sim_data", "p_$p,kpar_$kpar/")
    sim = Euler_integrator(system, 0.05, 1000, fps=30, Tplot=20, plot_functions=(plot_polymers!, plot_nematic_directors!, plot_velocity_vectors!), plotdim=2, Tsave=nothing, save_functions=(save_2d_polymer_polar_p!,), save_folder_path = save_folder); 
    return sim;

end 


#@sync @distributed

#=
for p in [0.04, 0.06, 0.08, 0.1, 0.13, 0.15, 0.2, 0.4]
    for (kpar, kperp) in [(-1., 0.), (1., 0.), (0., 1.), (0., -1.), (1/sqrt(2), 1/sqrt(2)),(-1/sqrt(2), 1/sqrt(2)),(1/sqrt(2), -1/sqrt(2)),(-1/sqrt(2), -1/sqrt(2))]
        display(p)
        display((kpar, kperp))
        sim = simulation(p, kpar, kperp)
    end
end
=#



sim = simulation(.2, -1., .5) 

# @profview sim = simulation() 

# @profview_allocs sim = simulation()  


 