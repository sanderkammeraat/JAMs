using Pkg
Pkg.add.(["StaticArrays","Distributions","Observables","GLMakie","JLD2","CodecZlib","HDF5","Distributed", "SlurmClusterManager","ProgressMeter"])

print(Threads.nthreads())

print("Now including Engine")
include(joinpath("..","src","Engine.jl"))

function relaxation_step(save_folder_path; Tsave=1000, Tplot=nothing)

    external_forces = []#[thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = (soft_disk_force(1,1.),)
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    N=1000
    ϕ = 1.
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.01], [0.001], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


   
    sizes =  [L, L, 4]


    initial_field_state=[]
    field_forces = []
    field_updaters = []

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, true,2.5*(1+poly));

    #Run integration
    sim = Euler_integrator(system,5e-2, 5e3,Tsave=Tsave, Tplot=Tplot, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path, save_tag="rx")#, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


function self_aligning_step(rx_step,J,v0, Dr, seed,save_folder_path; Tsave=1000, Tplot=nothing)

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,J),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = (soft_disk_force(1,1.),)

    local_dofevolvers = [overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1)]
    global_dofevolvers = []
    field_dofevolvers = []

    sizes = rx_step.system.sizes
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    initial_particle_state =deepcopy(rx_step.final_particle_state)

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)
        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0
        p_i.Dr[1] = Dr
        p_i.v0[1] = v0

    end
    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,true,rx_step.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,1e-2, 1e4,Tsave=Tsave,seed=seed, Tplot=Tplot, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path, save_tag="sa")# , fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end


function relax_again_step(sa_step, save_folder_path; Tsave=1000, Tplot=nothing)

    external_forces =[] # [thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = (soft_disk_force(1,1.),)

    local_dofevolvers = [overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1)]
    global_dofevolvers = []
    field_dofevolvers = []

    sizes = sa_step.system.sizes
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    initial_particle_state =deepcopy(sa_step.final_particle_state)

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)
        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0
        p_i.Dr[1] = 0.
        p_i.v0[1] = 0.

    end
    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,true,sa_step.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,5e-2, 5e3,Tsave=Tsave,seed=nothing, Tplot=Tplot, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path, save_tag="ra")# , fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end

#@profview rx_result= relaxation_step("",Tsave=nothing, Tplot=nothing)
#sa_result=self_aligning_step(rx_result,0.8, 0.0001, 0.00,1, ""; Tsave=nothing, Tplot=100);

#ra_result=relax_again_step(sa_result, ""; Tsave=nothing, Tplot=100);






# Drs = [0.,0.001, 0.01,0.02,0.05, 0.1, 0.2, 0.5, 1, 10] 
# Js=[0, 0.01, 0.1, 0.2, 0.5, 1. ,2., 5.]

Drs = [0.01] 
Js=[0.1]

seeds = reshape( collect(1:length(Drs)*length(Js)), (length(Drs),length(Js)) )

for j in eachindex(Js)
    for i in eachindex(Drs)

        J = Js[j]
        Dr = Drs[i]

        seed = seeds[i,j]
        print(seed)
        display("Running")

        base_path="/Volumes/T7_Shield/"
        #base_path = homedir()
        save_folder_path = joinpath(base_path,"sa","survey","free_disordered","phi_1","N_1000","t1e4","simdata", "J_$J","Dr_$Dr","seed_$seed");
        print(save_folder_path)

        rx_result = relaxation_step(save_folder_path)

        sa_result = self_aligning_step(rx_result,J, 0.0001,Dr, seed, save_folder_path);

        ra_result=relax_again_step(sa_result, save_folder_path);


    end
end
