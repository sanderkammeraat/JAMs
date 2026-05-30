#using Distributed
#using SlurmClusterManager
#addprocs(SlurmManager())

include(joinpath("..","src","Engine.jl"))


include(joinpath("..","io","InitialPositionGenerators.jl"))


function relaxation_step(save_folder_path; Tsave=nothing, Tplot=nothing,Nlin=20)

    external_forces = []#[thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = (soft_disk_force([1, 2],[1. 2.; 2. 1.]),)
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []
    #First make stair
    Nrows = 2*Nlin
    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]
    xs = []
    ys = []
    r=1.0
    ϕ=1.3
    l = 2

    #if regular hexagonal packing
    Nint_regular = (2 * (Nlin-1)-1)^2 - (Nlin-1)*((Nlin-1)-1)
    #corresponds to a packing fraction  pi/(2 sqrt(3))    
    #Total area available for interior particles is thus
    # Atot = Nint_regular * pi *  1^2/ (pi/(2 sqrt(3))) 
    Atot = Nint_regular * 1 * 2 * sqrt(3)

    Nint = ceil(Int64, ϕ * Atot/(pi*1^2) )


    typess = []
    #Setup a regular hexagonal lattice
    for row in 1:Nlin

        push!(xs, [xi-(Nlin+1)*l/2-(row-1)*l/2 for xi in l.*range(1,Nlin+row-1) ])

        push!(ys, [-(Nlin-row)/2*l*sqrt(3) for n in xs[row] ])

        if row==1
            rowtypes=2*ones(Int64, Nlin)

        else
            rowtypes = ones(Int64,length(xs[row]))
            rowtypes[1]=2
            rowtypes[end]=2

        end
        push!(typess,rowtypes )
    end

    for row in 1:Nlin-1
        push!(xs, xs[Nlin-row])
        push!(ys, ys[Nlin][1].-ys[Nlin-row])
        push!(typess, typess[Nlin-row])
    end

    xregular = vcat(xs...)
    yregular = vcat(ys...)
    types= vcat(typess...)
    Nregular = length(xregular)
    Nb = 6*(Nlin-1)

    N = Nint+Nb

    poly=0.15
    Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)
    while mean(Rs)<1 || mean(Rs)>1+ 1e-2
        Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)
        println(mean(Rs))
    end
    #we can input all coordinates, since they are only used to extract the min and max
    xi, yi = Random_in_hexagon(Nint, xregular, yregular, Rs)

    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]
    id=1
    #First loop over desired  normal particles

    for i=1:Nint
        push!(initial_state, PolarParticle3d([id], [1], [1], [1], [Rs[id]], [0.], [0.], [xi[i],yi[i],0.],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
        id+=1
    end

    #only pick out the boundary points of the regular hexagon
    for i=1:Nregular

        if types[i]==2
            push!(initial_state,ConfinedPolarParticle3d([id],[2], [1],[1], [Rs[id]], [0.], [0.], [xregular[i] , yregular[i],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]))
            id+=1
        end
    end

    size = [Nrows*l+2*l,Nrows*l+2*l,1];
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, false,2.5*r*(1+poly));

    #Run integration
    sim = Euler_integrator(system,0.05, 1e3,Tsave=Tsave, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path,save_tag="rx", fps=Inf, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tplot=Tplot); 
    return sim

end

function self_aligning_step(rx_step,J,v0, Dr, seed,save_folder_path; Tsave=nothing, Tplot=nothing)

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,J),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = (soft_disk_force([1, 2],[1 2; 2 1]),)

    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
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
    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,false,rx_step.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,0.01, 1e4,Tsave=Tsave,seed=seed, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path,save_tag="sa", fps=Inf, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tplot=Tplot); 
    return sim
end


function relax_again_step(sa_step, save_folder_path; Tsave=nothing, Tplot=nothing)

    external_forces =[] # [thermal_translational_noise(1, 0 .*[1.,1.,0])]

    pair_forces = (soft_disk_force([1, 2],[1. 2.; 2. 1.]),)

    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
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
    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers,field_dofevolvers,false,sa_step.system.rcut_pair_global);

    #Run integration
    sim = Euler_integrator(system,0.05, 1e3,Tsave=Tsave,seed=nothing, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path,save_tag="ra", fps=Inf, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2, Tplot=Tplot)# , fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim
end







# rx_result= relaxation_step(save_folder_path, Tsave=2000)

# # rx_step,J,v0, Dr, seed, save_folder_path
# sa_result=self_aligning_step(rx_result, .1 , 0.001, 0.01,100, save_folder_path, Tplot=100,  Tsave=nothing);

# ra_result=relax_again_step(sa_result, save_folder_path);

base_path = "/Volumes/T7_shield/sa/statistics"

Drs = [0.01, 0.1] 
Js=[0.1, 0.01]
Nlins = (30, 40, 50, 60)

#Let's keep the seeds constant across parameter choice
seeds = collect(1:10)

for k in eachindex(Nlins)

    Nlin=Nlins[k]
    #Creating the seeds
    for m in eachindex(seeds)
        seed = seeds[m]

        #We use the same 10 random inital conditions for all parameter space points, this to make sure the modes are the same

        display("Relaxation step")

        rx_save_folder_path = joinpath(base_path,"hex_disordered","phi_1.3","vary_Nlin","Nlin_$(Nlin)","simdata","initial_relaxation","seed_$seed")
        rx_result= relaxation_step(rx_save_folder_path, Tsave=1000,Nlin=Nlin)

        

        for i in eachindex(Drs)


                display(Threads.nthreads())
                println(Threads.nthreads())

                J = Js[i]
                Dr = Drs[i]
                v0 = 0.01
                seed = seeds[m]
                
                

                save_folder_path = joinpath(base_path,"hex_disordered","phi_1.3","vary_Nlin","Nlin_$(Nlin)","simdata","Dr_$Dr","J_$J","seed_$seed");
                display(save_folder_path)

                # rx_step,J,v0, Dr, seed, save_folder_path
                sa_result=self_aligning_step(deepcopy(rx_result), J , v0, Dr,nothing, save_folder_path,Tsave=100);

                ra_result=relax_again_step(deepcopy(sa_result), save_folder_path, Tsave=1000);

        end
    end
end

# for m in eachindex(seeds)
#     seed = seeds[m]

#     rx_save_folder_path = joinpath(base_path,"hex_disordered","phi_1.3","Nlin_20","simdata","initial_relaxation","seed_$seed","rx_JAMs_final_state.jld2")


#     jldopen(rx_save_folder_path, "r", iotype=IOStream ) do rx_file


#         rx_result = rx_file["SIM"]


#         for k in eachindex(v0s)

#             for i in eachindex(Drs)

#                     @sync @distributed for j in eachindex(Js)
                


#                     display(Threads.nthreads())

#                     J = Js[j]
#                     Dr = Drs[i]
#                     v0 = v0s[k]
#                     seed = seeds[m]
                    
                    

#                     save_folder_path = joinpath(base_path,"hex_disordered","phi_1.3","Nlin_20","simdata","v0_$(v0)","Dr_$Dr","J_$J","seed_$seed");
#                     display(save_folder_path)

#                     # rx_step,J,v0, Dr, seed, save_folder_path
#                     sa_result=self_aligning_step(deepcopy(rx_result), J , v0, Dr,nothing, save_folder_path,Tsave=100);

#                     ra_result=relax_again_step(deepcopy(sa_result), save_folder_path, Tsave=1000);


#                 end
#             end
#         end
#     end
# end