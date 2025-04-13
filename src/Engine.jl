using StaticArrays
using SparseArrays
using Observables
using JLD2
using CodecZlib
using ProgressMeter
#Note, only arrays can be changed in a struct. So initializing a struct attribute as array allows to change
#Type declaration in structs is important for performance, see https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations
include("Particles.jl")
include("Fields.jl")
include("Forces.jl")
include("DOFevolvers.jl")
include("FieldUpdaters.jl")
try # using try-catch to avoid error when there is not DISPLAY
    include("LivePlottingFunctions.jl")
catch LoadError
    @warn "Cannot load \"LivePlottingFunctions.jl\"."
end
include("SaveFunctions.jl")

@views function periodic!(p_i, systemsizes)

    for (i, xi) in pairs(p_i.x)

        if xi<-systemsizes[i]/2
            p_i.x[i] = xi + systemsizes[i]
        elseif xi>=systemsizes[i]/2
            p_i.x[i] = xi - systemsizes[i]
        end
    end
    return p_i
end

#Initialize unwrapped coordinates to save the user the hassle to set equal to the initial wrapped coordinates
function init_unwrap!(p_i, t)

    if t==0
        p_i.xuw.= copy(p_i.x)
    end

    return p_i
end


@views function minimal_image_closest_bin_center!(field_indices,x, bin_centers,system_sizes,system_Periodic)

    for (i, xi) in pairs(x)

        di = abs.( bin_centers[i].- x[i])
        field_indices[i] = argmin(di)

   
    end
    return field_indices

end
#dx is assumed to be pre-allocated
function minimal_image_difference!(dx,xi, xj, system_sizes, system_Periodic)

    
    for n in eachindex(xi)
        dx[n]=xj[n]-xi[n]
        if system_Periodic
            if dx[n]>system_sizes[n]/2
                dx[n]-=system_sizes[n]
            end
            if dx[n]<=-system_sizes[n]/2
                dx[n]+=system_sizes[n]
            end
        end 
    end
    return dx
end

"""
Description of the physical system in JAMs.

sizes: Array containing the linear size of the system. 

initial_particle_state: Array containing (struct instances of) particles

initial_field_state: Array containing (struct instances of) fields

external_forces: Array containing (struct instances of) forces that are calculated without information of the other particles or fields in the system.

pair_forces: Array containing (struct instances of) forces between two particles

field_forces: Array containing (struct instances of) forces acting on particles by a field

field_updaters: Array containing (struct instances of) forces/update rules acting on the field

dofevolvers: Array containing (struct instances of) evolvers that update the degrees of freedom of particles and fields

periodic: Bool describing whether the system is spatially periodic of the system sizes.

rcut_pair_global: Float setting the cutoff of all pair_forces and is used to generate cell lists.

"""
struct System{T1,T2, T3, T4, T5, T6, T7, T8, T9}

    #Vector that determines the linear size of the system
    sizes::T1

    #Array containing particles in a specific state
    initial_particle_state::T2

    #Array containing fields in a specific state
    initial_field_state::T3

    #Array of force functions:
    external_forces::T4

    pair_forces::T5

    field_forces::T6

    field_updaters::T7
    
    #Array of functions to evolve dof (and reinitialize forces)
    local_dofevolvers::T8

    global_dofevolvers::T9

    #Spatially periodic boundary conditions?
    Periodic::Bool

    #Global cutoff for pairwise interactions
    rcut_pair_global::Float64
end
#Output formatter to conveniently chain simulations in one .jl file without the need of intermediate saving to disk
struct SIM{T1, T2, T3, T4}
    final_particle_state::T1
    final_field_state::T2
    dt::T3
    t_stop::T4
    system::System
end
function save_raw_force_data!(file, preamble, force)

    force_name = string(nameof(typeof(force)))

    field_names = fieldnames(typeof(force))

    for field_name in field_names

        name = string(field_name)
        val = getfield(force, field_name)


        file[preamble*force_name*"/"*name] = val

    end

    return file
end

function save_raw_metadata!(file, system, integration_tax,dt,t_stop,Tsave,save_tax, master_seed)

    for force in system.external_forces

        preamble = "system/forces/external_forces/"

        save_raw_force_data!(file,preamble, force)

    end
    for force in system.pair_forces

        preamble = "system/forces/pair_forces/"

        save_raw_force_data!(file,preamble, force)

    end
    for force in system.field_forces

        preamble = "system/forces/field_forces/"

        save_raw_force_data!(file,preamble, force)

    end
    for fieldupdater in system.field_updaters

        preamble = "system/field_updaters/"

        save_raw_force_data!(file,preamble, fieldupdater)

    end

    for dofevolver in system.local_dofevolvers

        dofevolver_name = string(nameof(dofevolver))
        file["system/local_dofevolvers/"*dofevolver_name] = dofevolver_name
    end

    for dofevolver in system.global_dofevolvers

        dofevolver_name = string(nameof(dofevolver))
        file["system/global_dofevolvers/"*dofevolver_name] = dofevolver_name
    end

    file["integration_info/integration_tax"] = integration_tax

    file["integration_info/Tsave"] = Tsave

    file["integration_info/save_tax"] = save_tax

    file["integration_info/dt"] = dt

    file["integration_info/t_stop"] = t_stop

    file["integration_info/master_seed"] = master_seed

    return file

end



function Euler_integrator(system, dt, t_stop; seed=nothing, Tsave=nothing, save_functions=nothing, save_folder_path=nothing, save_tag=nothing, Tplot=nothing, fps=nothing, plot_functions=nothing,plotdim=nothing)


    integration_tax = collect(0:dt:t_stop)



    #This is really amazing: Julia rng is dependent on the task spawn structure, NOT on the 
    # parallel executation schedule, see https://julialang.org/blog/2021/11/julia-1.7-highlights/#new_rng_reproducible_rng_in_tasks
    # for more info, so this makes the program, even with mulitithreading and dynamic thread scheduling reproducible if we create a rng for every particle and field!
    # See for info https://github.com/JuliaLang/julia/issues/49064 !!
    if !isnothing(seed)
        master_seed = seed
        Random.seed!(seed)

    else
    #Otherwise generate one for reproducibility
    #Reset seed in case user used a specific seed for initial conditions
        Random.seed!()
        seed = rand(1:100000000000000000000000000000)
        master_seed = seed
        Random.seed!(seed)
    end
    if !isnothing(Tsave)
        save_nax = [n for n in eachindex(integration_tax) if (n-1)%Tsave==0]
        save_tax = [ integration_tax[n] for n in eachindex(integration_tax) if (n-1)%Tsave==0 ]

        #Prepare save folder
        #If folder already exists, simply returns folder path. If folder is not existing, it will create the (sub)folders and return the path
        mkpath(save_folder_path)

        if isnothing(save_tag)
            JAMs_file_name = "JAMs_container.jld2"
            
            raw_data_file_name = "raw_data.jld2"
        else
            JAMs_file_name = save_tag * "_"* "JAMs_container.jld2"
            
            raw_data_file_name = save_tag * "_"*"raw_data.jld2"
        end

        if isfile(joinpath(save_folder_path, JAMs_file_name)) || isfile(joinpath(save_folder_path, raw_data_file_name))
            error("JAMs: Specified save folder already contains JAMs file(s): " * JAMs_file_name * " and/or " * raw_data_file_name*". JAMs aborted to prevent overwriting.")
        end

        #Store system and integration info in JAMS container
        jldopen(joinpath(save_folder_path, JAMs_file_name),"a+") do JAMs_file

            JAMs_file["system"] = system

            JAMs_file["integration_info/integration_tax"] = integration_tax

            JAMs_file["integration_info/dt"] = dt

            JAMs_file["integration_info/t_stop"] = t_stop

            JAMs_file["integration_info/Tsave"] = Tsave

            JAMs_file["integration_info/save_tax"] = save_tax

            JAMs_file["integration_info/master_seed"] = master_seed

            if !isnothing(save_functions)
                JAMs_file["integration_info/save_functions"] = save_functions
            end

            if !isnothing(Tplot)

                JAMs_file["integration_info/Tplot"] = Tplot

                JAMs_file["integration_info/fps"] = fps

                JAMs_file["integration_info/plot_functions"] = plot_functions

                JAMs_file["integration_info/plotdim"] = plotdim


            end
            


        end

        #Store similar info in raw_data container
        jldopen(joinpath(save_folder_path, raw_data_file_name),"a+") do raw_data_file

            save_raw_metadata!(raw_data_file, system, integration_tax,dt, t_stop, Tsave,save_tax,master_seed)

        end
   
    end

    if !isnothing(Tsave)
        n_final_save = save_nax[end]
    else
        n_final_save = length(integration_tax)
    end

    if system.Periodic==false
        @warn ("JAMs: System is set to non-periodic: you should make sure particles always stay in system sizes for correctly working cell lists. The program will catch this by throwing an error if a particle is detected outside the box.")
    end

    system, cells,cell_bin_centers,stencils = construct_cell_lists!(system)

    Npair = length(system.pair_forces)

    Nfield = length(system.field_forces)

    Nfieldu = length(system.field_updaters)

    

    current_particle_state = copy(system.initial_particle_state)
    current_field_state = copy(system.initial_field_state)

    #We do not really care that  it is a shallow copy, we will properly set it at the end. We need to prepare it as variables that exist outside the for loop over time.
    final_particle_state = copy(system.initial_particle_state)
    final_field_state = copy(system.initial_field_state)

    rngs_fields = [Xoshiro(master_seed+i) for i in eachindex(current_field_state)]

    #Assuming number of fields stay constant
    rngs_particles = [Xoshiro(length(current_field_state)+master_seed+i) for i in eachindex(current_particle_state)]

    if !isnothing(Tplot)
        if fps!=0
            cpsO = Observable(current_particle_state)
            cfsO = Observable(current_field_state)
            tO = Observable(0.)
            f, ax = setup_system_plotting(system.sizes,plot_functions, plotdim,cpsO,cfsO,tO)
        end
    end


    #Loop over time
    #Variable to keep track of the number of frames saved
    frame_counter = 1
    @showprogress dt = 1 desc="JAMming in progress..." showspeed=true for (n, t) in pairs(integration_tax)
        #Threads.@threads 
        #Loop over particles and write generalized forces to p_i in place!
        Threads.@threads for i in eachindex(current_particle_state)
            p_i = current_particle_state[i]
            p_i=init_unwrap!(p_i, t)
            p_i = particle_step!(i,p_i, current_particle_state,Npair,t, dt, system,cells,cell_bin_centers,stencils,rngs_particles)
            current_particle_state[i]=p_i
        end


        #Threadsafe field forces
        if Nfield>0
            for i in eachindex(current_particle_state)
                p_i = current_particle_state[i]
                p_i, current_field_state = contribute_field_forces!(p_i, current_field_state, t, dt,system,rngs_particles)
                current_particle_state[i]=p_i
            end
        end

        #Field update loop
        for i in eachindex(current_field_state)

            field_i = current_field_state[i]

            field_i = field_step!(i, field_i,current_field_state,Nfieldu,t,dt, system, rngs_fields)

            current_field_state[i] = field_i
        end

        if !isnothing(Tsave) && !isnothing(save_functions)
            if (n-1)%Tsave==0

                jldopen(joinpath(save_folder_path, raw_data_file_name),"a+") do raw_data_file
 
                    for save_function in save_functions

                        raw_data_file=save_function(raw_data_file,current_particle_state,current_field_state, n, Tsave, t,frame_counter)
                        
                    end
                end

                frame_counter+=1
            end
        end
        #Save the states before the final dof step

        if n==n_final_save
            if !isnothing(Tsave)
                jldopen(joinpath(save_folder_path, JAMs_file_name),"a+") do JAMs_file

                    JAMs_file["final_particle_state"] = current_particle_state
        
                    JAMs_file["final_field_state"] = current_field_state
        
                end
            end
            final_particle_state = deepcopy(current_particle_state)
            final_field_state = deepcopy(current_field_state)
        end


        #Only now evolve dofs of every particle
        #Local
        Threads.@threads for i in eachindex(current_particle_state)
            p_i = current_particle_state[i]
            for dofevolver in system.local_dofevolvers
                p_i=evolve_locally!(p_i, t, dt, dofevolver)
            end
        end
        #Or global

        for dofevolver in system.global_dofevolvers
            current_particle_state,current_field_state =evolve_globally!(current_particle_state, current_field_state, system, cells, stencils, dt, dofevolver)
        end

        #Perform checks
        Threads.@threads for i in eachindex(current_particle_state)
            p_i = current_particle_state[i]

            #apply periodic boundary conditions
            if system.Periodic
                p_i=periodic!(p_i, system.sizes)
        
            else
                check_outside_system(p_i,system.sizes)
            end
            current_particle_state[i] = p_i
        end


        current_particle_state,cells = update_cells!(current_particle_state, cells, cell_bin_centers,system)

        #DOF evolver fields
        for i in eachindex(current_field_state)

            field_i = current_field_state[i]

            for dofevolver in system.dofevolvers
                field_i=dofevolver(field_i, t, dt)
            end
        end


        cells = update_ghost_cells!(cells,system)

        
        if !isnothing(Tplot)
            if fps!=0
                if (n-1)%Tplot==0
                    cpsO[] = current_particle_state
                    cfsO[]= current_field_state
                    tO[] = t
                    sleep(1/fps)
                end
            end
        end

    end
    return SIM(final_particle_state, final_field_state, dt, t_stop, system);
end




function particle_step!(i,p_i, current_particle_state,Npair,t, dt, system,cells,cell_bin_centers,stencils,rngs_particles)
    if Npair>0
        p_i=contribute_pair_forces!(i,p_i, current_particle_state,t, dt, system,cells,stencils,rngs_particles)
    end

    for force in system.external_forces
        p_i=contribute_external_force!(p_i, t, dt, force,rngs_particles)
    end

    return p_i
end

function check_outside_system(p_i, system_sizes)
    for i in eachindex(p_i.x)
        if p_i.x[i]>system_sizes[i]/2 || p_i.x[i]<-system_sizes[i]/2
            error("JAMs: Particle outside simulation box. This invalidates cell lists. Suggested fix: make sure particles always stay inside system sizes by increasing the system size of the relevant dimension.")
        end
    end

end

#currently not using current_field_state as there are currently no pair interactions betwen fields, but allow for the implementation.
function field_step!(i, field_i,current_field_state,Nfieldu,t,dt, system, rngs_fields)
    if Nfieldu>0
        for field_updater in system.field_updaters
            field_i=contribute_field_update!(field_i, t, dt, field_updater, rngs_fields)
        end
    end

    return field_i

end


function contribute_field_forces!(p_i, current_field_state, t, dt,system, rngs_particles)
    field_indices = @MVector zeros(Int,length(p_i.x))
    

    for force in system.field_forces

        for (j, field_j) in pairs(current_field_state)
            
            field_indices = minimal_image_closest_bin_center!(field_indices, p_i.x, field_j.bin_centers, system.sizes, system.Periodic)

            p_i, current_field_state[j] = contribute_field_force!(p_i, field_j, field_indices, t, dt, force,rngs_particles)
        end
    end
    return p_i, current_field_state
end




function contribute_pair_forces!(i,p_i, current_particle_state, t, dt,system,cells,stencils,rngs_particles)
    
    dx = @MVector zeros(Float64,length(p_i.x))
    neighbours= get_neighbours(p_i,cells,stencils)
    if !isnothing(neighbours)

        for n in neighbours

            if i!=n
                p_j = current_particle_state[n]

                dx = minimal_image_difference!(dx, p_i.x, p_j.x, system.sizes, system.Periodic)

                dxn = norm(dx)
                
                if dxn<=system.rcut_pair_global
                    for force in system.pair_forces
                        p_i=contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force,rngs_particles)
                    end
                end

            end
        end
    end
    return p_i
end



function construct_cell_lists!(system)

    dim = length(system.sizes)
    lbin = system.rcut_pair_global

    Lx = system.sizes[1]
    Ly = system.sizes[2]
    #First and last bin center should be far outside the simulation box so a particle is never associated in a ghost cell
    #via the minimal_image_closest_bin_center function. Note that system spans from -Li/2 to Li/2 so -Li should be
    #safely outside reach both when using periodic boundary conditions
    #or when using a finite system, because particle outside box will raise and error and  abort program
    #In case system size is entered as Int64, without FLoat64 in front, the array will be typed as Int64,
    # not allowing appending the  float values in the next line. Therefore declare type as Float64[].
    # +lbin incase lbin is larger than one of the system sizes
    x_bin_centers = Float64[-Lx-lbin]
    x_bin_centers = append!(x_bin_centers,range(start=-Lx/2, stop=Lx/2, step=lbin).+lbin/2)
    x_bin_centers = append!(x_bin_centers,Lx+lbin)

    y_bin_centers = Float64[-Ly-lbin]
    y_bin_centers = append!(y_bin_centers,range(start=-Ly/2, stop=Ly/2, step=lbin).+lbin/2)
    y_bin_centers = append!(y_bin_centers,Ly+lbin)
    nx = length(x_bin_centers)
    ny = length(y_bin_centers)
    if dim==2
        cell_bin_centers = [x_bin_centers, y_bin_centers]

        cells = reshape([Int64[] for i in 1:nx*ny],nx,ny)
        for p_i in system.initial_particle_state

            cell_indices=@MVector zeros(Int,length(p_i.x))
            cell_indices = minimal_image_closest_bin_center!(cell_indices,p_i.x, cell_bin_centers,system.sizes,system.Periodic)

            #push!(Int64[id for id in  cells[cell_indices][1]]))
            cells[cell_indices...]= append!(cells[cell_indices...],p_i.id[1])
            p_i.ci.= cell_indices

        end
        stencils = [ @SVector [ni, nj] for ni in -1:1 for nj in -1:1]
        

    elseif dim==3
        Lz = system.sizes[3]
        z_bin_centers = Float64[-Lz-lbin]
        z_bin_centers = append!(z_bin_centers,range(start=-Lz/2, stop=Lz/2, step=lbin).+lbin/2)
        z_bin_centers = append!(z_bin_centers,Lz+lbin)
        nz = length(z_bin_centers)
        cell_bin_centers = [x_bin_centers, y_bin_centers, z_bin_centers]

        cells = reshape([Int64[] for i in 1:nx*ny*nz],nx,ny,nz)
        for p_i in system.initial_particle_state

            cell_indices=@MVector zeros(Int,length(p_i.x))
            cell_indices = minimal_image_closest_bin_center!(cell_indices,p_i.x, cell_bin_centers,system.sizes,system.Periodic)
            

            #push!(Int64[id for id in  cells[cell_indices][1]]))
            cells[cell_indices...]= append!(cells[cell_indices...],p_i.id[1])
            p_i.ci.= cell_indices

        end
        stencils = [ @SVector [ni, nj, nk] for ni in -1:1 for nj in -1:1 for nk in -1:1]
    end
    @views cells = update_ghost_cells!(cells,system)

    return system, cells, cell_bin_centers, stencils
end

function get_neighbours(p_i, cells, stencils)
    candidate_cell_ind=@MVector zeros(Int64, length(p_i.ci))
    n=0
    #First check number
    for stencil in stencils

        for i in eachindex(candidate_cell_ind)
            candidate_cell_ind[i]= p_i.ci[i] + stencil[i]
        end

        if !isempty(cells[candidate_cell_ind...])
            for id in cells[candidate_cell_ind...]
                #push!(neighbours,id)
                n+=1
            end
        end
    end
    #Then allocate
    if n>0
        neighbours = zeros(Int64, n)
        m=1
        for stencil in stencils

            for i in eachindex(candidate_cell_ind)
                candidate_cell_ind[i]= p_i.ci[i] + stencil[i]
            end
    
            if !isempty(cells[candidate_cell_ind...])
                for id in cells[candidate_cell_ind...]
                    neighbours[m]=id
                    m+=1
                end
            end
        end
    else
        neighbours = nothing
    end

    return neighbours
end


function find_new_bin_locations(current_particle_state, cell_bin_centers,system)

    lbin = system.rcut_pair_global

    #Collect old locations to preallocate for new one
    new_bin_locations = [copy(p_i.ci) for p_i in current_particle_state]
    Threads.@threads for i in eachindex(current_particle_state)

        p_i = current_particle_state[i]

        for j in eachindex(p_i.ci)

            new_bin_locations[i][j]+= round(Int64, (p_i.x[j] - cell_bin_centers[j][p_i.ci[j]])/lbin )

        end

    end

    return new_bin_locations

end

function update_cells!(current_particle_state, cells, cell_bin_centers,system)

    new_bin_locations=find_new_bin_locations(current_particle_state, cell_bin_centers,system)

    #Updating the cell list is not threadsafe, hence put it outside the threaded loop
    #The cell list update must come *after* the DOF evolving of particles!
    for i in eachindex(current_particle_state)

        p_i = current_particle_state[i]
        new_bin_location = new_bin_locations[i]
        if p_i.id[1] in cells[new_bin_location...]
        else
            #remove
            filter!(e->e≠p_i.id[1],cells[p_i.ci...])
            #add to correct lists
            p_i.ci.= new_bin_location

            push!(cells[new_bin_location...],p_i.id[1])
        end
        current_particle_state[i] = p_i
    end
    return current_particle_state,cells
end



@views function update_ghost_cells!(cells,system)
    
    if system.Periodic
        dims = length(system.sizes)
        if dims==2
            cells[1,:].= cells[end-1,:]
            cells[end,:].= cells[2,:]

            cells[:,1].= cells[:,end-1]
            cells[:,end].= cells[:,2]

            cells[1,1] .= cells[end-1, end-1]
            cells[1,end] .= cells[end-1, 2]

            cells[end,1] .= cells[2, end-1]
            cells[end,end] .= cells[2, 2]

        end

        if dims==3

            cells[1,:,:].= cells[end-1,:,:]
            cells[end,:,:].= cells[2,:,:]

            cells[:,1,:].=  cells[:,end-1,:]
            cells[:,end,:].= cells[:,2,:]

            cells[:,:,1].=  cells[:,:,end-1]
            cells[:,:,end].=  cells[:,:,2]

            #set the rings
            #Write for three faces sharing a vertex, comment out duplicates, then do the remaining edges
            ##face 1
            cells[1,:,1].= cells[end-1,:,end-1]
            cells[1,:,end].= cells[end-1,:,2]

            cells[1,1,:].= cells[end-1,end-1,:]
            cells[1,end,:].= cells[end-1,2,:]

            ## face 2

            cells[:,1,1].= cells[:,end-1, end-1]
            cells[:,end,1].= cells[:,2,end-1]

            #cells[1,:,1] = cells[end-1,:,end-1]
            cells[end,:,1].= cells[2,:,end-1]

            ## face 3

            #cells[1,1,:] = cells[end-1, end-1,:]
            cells[end,1,:].= cells[2,end-1,:]

            #cells[:,1,1] = cells[:,end-1,end-1]
            cells[:,1,end].= cells[:,end-1,2]

            #We are now at 9 of 12 edges
            cells[end, end,:].= cells[2,2,:]
            cells[end,:,end].= cells[2,:,2]
            cells[:,end, end].= cells[:,2,2]
            #Done

            #set the 8 corner points

            cells[1,1,1].=cells[end-1, end-1, end-1]

            cells[end, 1, 1] .=cells[2, end-1, end-1]
            cells[1,1,end] .=cells[end-1, end-1, 2]
            cells[1,end,1] .= cells[end-1,2,end-1]

            cells[1, end, end] .= cells[end-1, 2, 2]
            cells[end, end, 1] .= cells[2,2, end-1]
            cells[end, 1, end] .= cells[2, end-1, 2]

            cells[end, end, end] .= cells[2,2,2]

        end

    end
    return cells
end
