include(joinpath("../src","Engine.jl"))

include(joinpath("../io","InitialPositionGenerators.jl"))

# @everywhere include(joinpath("..","src","Engine.jl"))


# @everywhere include(joinpath("..","io","InitialPositionGenerators.jl"))


function relaxation_step(save_folder_path; Tsave=nothing, Tplot=nothing, N=2000)

    external_forces = []#[thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = (soft_disk_force(1,1.),)
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []
    N=N
    ϕ = 1.
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  sqrt(pi *sum(Rs.^2) / ϕ)

    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.01], [0.001], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


    sizes = [L,L,1];
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, true,2.5*(1+poly));

    #Run integration
    sim = Euler_integrator(system,0.05, 1e3,Tsave=Tsave, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path,save_tag="rx")#, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tplot=Tplot); 
    return sim

end



function self_aligning_step(rx_step,J,v0, Dr, seed,save_folder_path; Tsave=nothing, Tplot=nothing)

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
    sim = Euler_integrator(system,0.01, 1e4,Tsave=Tsave,seed=seed, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path,save_tag="sa")#,  plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tplot=Tplot); 
    return sim
end


function relax_again_step(sa_step, save_folder_path; Tsave=nothing, Tplot=nothing)

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
    sim = Euler_integrator(system,0.05, 1e3,Tsave=Tsave,seed=nothing, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path,save_tag="ra")#, fps=Inf, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tplot=Tplot)# , fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end







# rx_result= relaxation_step(save_folder_path, Tsave=2000)

# # rx_step,J,v0, Dr, seed, save_folder_path
# sa_result=self_aligning_step(rx_result, .1 , 0.001, 0.01,100, save_folder_path, Tplot=100,  Tsave=nothing);

# ra_result=relax_again_step(sa_result, save_folder_path);

#@everywhere base_path = "/Users/kammeraat/test_free/"#"/data1/kammeraatsc1/sa/statistics/free/"


base_path ="/Volumes/T7_Shield/sa/statistics/free/"


Drs = [ 0.01, 0.1] 

#13 Js
Js=[0.0, 0.001, 0.01, 0.02, 0.04, 0.08, 0.10, 0.12, 0.16, 0.2, 0.4, 0.8, 1]#[0.1, 0.01]
v0s = [0.01]#[0.001, 0.003, 0.01, 0.03,  0.10, 0.12, 0.15]

#Let's keep the seeds constant across parameter choice
seeds = collect(1:10)

#Creating the seeds
for m in eachindex(seeds)
    seed = seeds[m]

    #We use the same 10 random inital conditions for all parameter space points, this to make sure the modes are the same

    display("Relaxation step")

    rx_save_folder_path = joinpath(base_path,"phi_1.0","N_2000","simdata","initial_relaxation","seed_$seed")
    rx_result= relaxation_step(rx_save_folder_path, Tsave=1000)

    for k in eachindex(v0s)

        for i in eachindex(Drs)

                for j in eachindex(Js)
             

                display(Threads.nthreads())

                J = Js[j]
                Dr = Drs[i]
                v0 = v0s[k]
                seed = seeds[m]
                
                

                save_folder_path = joinpath(base_path,"phi_1.0","N_2000","simdata","v0_$(v0)","Dr_$Dr","J_$J","seed_$seed");
                display(save_folder_path)

                # rx_step,J,v0, Dr, seed, save_folder_path
                sa_result=self_aligning_step(deepcopy(rx_result), J , v0, Dr,nothing, save_folder_path,Tsave=200);

                ra_result=relax_again_step(deepcopy(sa_result), save_folder_path, Tsave=1000);


            end
        end
    end
end

# for m in eachindex(seeds)
#     seed = seeds[m]

#     rx_save_folder_path = joinpath(base_path,"hex_disordered","phi_1.3","Nlin_20","simdata","initial_relaxation","seed_$seed","rx_JAMs_final_state.jld2")


#     jldopen(rx_save_folder_path, "r", iotype=IOStream ) do rx_file


#         rx_result = rx_file["SIM"]

#             for i in eachindex(Drs)

#                 @sync @distributed for k in eachindex(v0s)


#                     display(Threads.nthreads())

#                     J = Js[i]
#                     Dr = Drs[i]

#                     v0 = v0s[k]
#                     seed = seeds[m]
                    
                    

#                     save_folder_path = joinpath(base_path,"hex_disordered","phi_1.3","Nlin_20","simdata","v0_$(v0)","Dr_$Dr","J_$J","seed_$seed");
#                     display(save_folder_path)

#                     # rx_step,J,v0, Dr, seed, save_folder_path
#                     sa_result=self_aligning_step(deepcopy(rx_result), J , v0, Dr,nothing, save_folder_path,Tsave=100);

#                     ra_result=relax_again_step(deepcopy(sa_result), save_folder_path, Tsave=1000);


#             end
#         end
#     end
# end