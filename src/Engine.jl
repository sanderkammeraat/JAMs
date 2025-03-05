using StaticArrays
using Observables
using JLD2
using CodecZlib
using ProgressMeter
using DataFrames
#Note, only arrays can be changed in a struct. So initializing a struct attribute as array allows to change
#Type declaration in structs is important for performance, see https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations
include("Particles.jl")
include("Fields.jl")
include("Forces.jl")
include("DOFevolvers.jl")
include("FieldUpdaters.jl")
include("LivePlottingFunctions.jl")
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



@views function minimal_image_closest_bin_center!(field_indices,x, bin_centers,system_sizes,system_Periodic)

    d = @MVector zeros(length(x))
    for (i, xi) in pairs(x)
        
        d= abs.( bin_centers[i].- x[i])
        field_indices[i] = argmin(d)
        
    end
    return field_indices

end

@views function minimal_image_difference!(dx,xi, xj, system_sizes, system_Periodic)

    
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
struct System{T0,T1, T2, T3, T4, T5, T6, T7}

    #Vector that determines the linear size of the system
    sizes::T0

    #Array containing particles in a specific state
    initial_particle_state::T1

    #Array containing fields in a specific state
    initial_field_state::T2

    #Array of force functions:
    external_forces::T3

    pair_forces::T4

    field_forces::T5

    field_updaters::T6
    
    #Array of functions to evolve dof (and reinitialize forces)
    dofevolvers::T7

    #Spatially periodic boundary conditions?
    Periodic::Bool

    #Global cutoff for pairwise interactions
    rcut_pair_global::Float64

end

struct SIM{T1, T2, T3, T4}
    final_particle_state::T1
    final_field_state::T2
    dt::T3
    t_stop::T4
    system::System
end


function Euler_integrator(system, dt, t_stop;  Tsave=nothing, save_functions=nothing, save_folder_path=nothing, Tplot=nothing, fps=nothing, plot_functions=nothing,plotdim=nothing)


    integration_tax = 0:dt:t_stop

    if !isnothing(Tsave)
        save_tax = [ integration_tax[n] for n in eachindex(integration_tax) if n%Tsave==0   ]

        #Prepare save folder
        mkpath(save_folder_path)

        JAMs_file_name = "JAMs_container.jld2"
        
        raw_data_file_name = "raw_data.jld2"


        #Store system in JAMS container
        jldopen(save_folder_path*JAMs_file_name,"a+") do JAMs_file

            JAMs_file["system"] = system

            JAMs_file["integration_tax"] = integration_tax

            JAMs_file["Tsave"] = Tsave

            JAMs_file["save_tax"] = save_tax

        end
    end

    if system.Periodic==false
        print("System is set to non-periodic: you should make sure particles always stay in system sizes for correctly working cell lists. The program will catch this by throwing an error if a particle is detected outside the box.")
    end

    system, cells,cell_bin_centers,stencils = construct_cell_lists!(system)

    Npair = length(system.pair_forces)

    Nfield = length(system.field_forces)

    Nfieldu = length(system.field_updaters)

    

    current_particle_state = copy(system.initial_particle_state)

    current_field_state = copy(system.initial_field_state)

    if !isnothing(plot_functions)
        if fps!=0
            cpsO = Observable(current_particle_state)
            cfsO = Observable(current_field_state)
            tO = Observable(0.)
            f, ax = setup_system_plotting(system.sizes,plot_functions, plotdim,cpsO,cfsO,tO)
        end
    end

    #Loop over time
    @showprogress dt = 1 desc="JAMming in progress..." showspeed=true for (n, t) in pairs(integration_tax)
        #Threads.@threads 
        #Loop over particles and write generalized forces to p_i in place!
        Threads.@threads for i in eachindex(current_particle_state)
            p_i = current_particle_state[i]
            p_i, cells,current_field_state = particle_step!(i,p_i, current_particle_state,current_field_state,Npair, Nfield,t, dt, system,cells,cell_bin_centers,stencils)
            current_particle_state[i]=p_i
        end

        #Field loop
        for i in eachindex(current_field_state)

            field_i = current_field_state[i]

            field_i = field_step!(i, field_i,current_field_state,Nfieldu,t,dt, system)

            current_field_state[i] = field_i
        end

        if !isnothing(Tsave) && !isnothing(save_functions)
            if n%Tsave==0

                jldopen(save_folder_path*raw_data_file_name,"a+") do raw_data_file

            
                    for save_function in save_functions

                        save_function(raw_data_file,current_particle_state,current_field_state, n, Tsave, t)

                    end

                end
            end
        end

        #Only now evolve dofs of every particle
        Threads.@threads for i in eachindex(current_particle_state)
            p_i = current_particle_state[i]
            for dofevolver in system.dofevolvers
                p_i=dofevolver(p_i, t, dt)
            end
            if system.Periodic
                p_i=periodic!(p_i, system.sizes)
        
            else
                check_outside_system(p_i,system.sizes)
            end
            current_particle_state[i] = p_i
        end

        #Updating the cell list is not threadsafe, hence put it outside the threaded loop
        for i in eachindex(current_particle_state)

            p_i = current_particle_state[i]
            p_i, cells=update_cells!(p_i, cells, cell_bin_centers, stencils,system)
            current_particle_state[i]=p_i
        end


        #DOF evolver fields
        for i in eachindex(current_field_state)

            field_i = current_field_state[i]

            for dofevolver in system.dofevolvers
                field_i=dofevolver(field_i, t, dt)
            end
        end


        cells = update_ghost_cells!(cells,system)

        
        if !isnothing(plot_functions)
            if fps!=0
                if n%Tplot==0
                    cpsO[] = current_particle_state
                    cfsO[]= current_field_state
                    tO[] = t
                    sleep(1/fps)
                end
            end
        end

    end

    if !isnothing(Tsave)
        jldopen(save_folder_path*JAMs_file_name,"a+") do JAMs_file

            JAMs_file["final_particle_state"] = current_particle_state

            JAMs_file["final_field_state"] = current_field_state

        end
    end

    return SIM(current_particle_state, current_field_state,dt, t_stop, system);
end
function particle_step!(i,p_i, current_particle_state,current_field_state,Npair, Nfield,t, dt, system,cells,cell_bin_centers,stencils)
    if Npair>0
        p_i=contribute_pair_forces!(i,p_i, current_particle_state,t, dt, system,cells,stencils)
    end

    if Nfield>0
        p_i, current_field_state = contribute_field_forces!(p_i, current_field_state, t, dt,system)

    end

    for force in system.external_forces
        p_i=contribute_external_force!(p_i, t, dt, force)
    end

    return p_i,cells, current_field_state
end

function check_outside_system(p_i, system_sizes)
    for i in eachindex(p_i.x)
        if p_i.x[i]>system_sizes[i]/2 || p_i.x[i]<-system_sizes[i]/2
            error("Particle outside simulation box. This invalidates cell lists. Suggested fix: make sure particles always stay inside system sizes by increasing the system size of the relevant dimension.")
        end
    end

end

#currently not using current_field_state as there are currently no pair interactions betwen fields, but allow for the implementation.
function field_step!(i, field_i,current_field_state,Nfieldu,t,dt, system)
    if Nfieldu>0
        for field_updater in system.field_updaters
            field_i=contribute_field_update!(field_i, t, dt, field_updater)
        end
    end

    return field_i

end

# function save_state!(states, current_state, n, Tsave)

#     if n%Tsave==0 && n>0
#         states=push!(states,deepcopy(current_state))

#     end
#     return states

# end

function contribute_field_forces!(p_i, current_field_state, t, dt,system)
    field_indices = @MVector zeros(Int,length(p_i.x))
    

    for force in system.field_forces

        for (j, field_j) in pairs(current_field_state)
            
            field_indices = minimal_image_closest_bin_center!(field_indices, p_i.x, field_j.bin_centers, system.sizes, system.Periodic)

            p_i, current_field_state[j] = contribute_field_force!(p_i, field_j, field_indices, t, dt, force)
        end
    end
    return p_i, current_field_state
end




function contribute_pair_forces!(i,p_i, current_particle_state, t, dt,system,cells,stencils)
    
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
                        p_i=contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force)
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
    x_bin_centers = [-1e8]
    x_bin_centers = append!(x_bin_centers,range(start=-Lx/2, stop=Lx/2, step=lbin).+lbin/2)
    x_bin_centers = append!(x_bin_centers,1e8)

    y_bin_centers = [-1e8]
    y_bin_centers = append!(y_bin_centers,range(start=-Ly/2, stop=Ly/2, step=lbin).+lbin/2)
    y_bin_centers = append!(y_bin_centers,1e8)
    nx = length(x_bin_centers)
    ny = length(y_bin_centers)
    if dim==2
        cell_bin_centers = [x_bin_centers, y_bin_centers]

        cells = reshape([Int64[] for i in 1:nx*ny],nx,ny)
        for p_i in system.initial_particle_state

            cell_indices=@MVector zeros(Int,length(p_i.x))
            cell_indices = minimal_image_closest_bin_center!(cell_indices,p_i.x, cell_bin_centers,system.sizes,system.Periodic)

            #push!(Int64[id for id in  cells[cell_indices][1]]))
            cells[cell_indices...]= append!(cells[cell_indices...],p_i.id)
            p_i.ci.= cell_indices

        end
        stencils = [ @SVector [ni, nj] for ni in -1:1 for nj in -1:1]
        

    elseif dim==3
        Lz = system.sizes[3]
        z_bin_centers = [-1e8]
        z_bin_centers = append!(z_bin_centers,range(start=-Lz/2, stop=Lz/2, step=lbin).+lbin/2)
        z_bin_centers = append!(z_bin_centers,1e8)
        nz = length(z_bin_centers)
        cell_bin_centers = [x_bin_centers, y_bin_centers, z_bin_centers]

        cells = reshape([Int64[] for i in 1:nx*ny*nz],nx,ny,nz)
        for p_i in system.initial_particle_state

            cell_indices=@MVector zeros(Int,length(p_i.x))
            cell_indices = minimal_image_closest_bin_center!(cell_indices,p_i.x, cell_bin_centers,system.sizes,system.Periodic)

            #push!(Int64[id for id in  cells[cell_indices][1]]))
            cells[cell_indices...]= append!(cells[cell_indices...],p_i.id)
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

@views function update_cells!(p_i, cells, cell_bin_centers, stencils,system)
    new_bin_location=@MVector zeros(Int,length(p_i.x))
    new_bin_location= minimal_image_closest_bin_center!(new_bin_location,p_i.x, cell_bin_centers,system.sizes,system.Periodic)
    if p_i.id in cells[new_bin_location...]
    else
        #find which ghost cell direction needs to be changed

        #remove
        filter!(e->e≠p_i.id,cells[p_i.ci...])
        #add to correct lists
        p_i.ci.= new_bin_location

        push!(cells[new_bin_location...],p_i.id)
    end
    return p_i,cells
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
