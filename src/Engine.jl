
using ProgressBars
using StaticArrays
using Observables
using JLD2
using CodecZlib

#Note, only arrays can be changed in a struct. So initializing a struct attribute as array allows to change
#Type declaration in structs is important for performance, see https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations
include("Particles.jl")
include("Fields.jl")
include("Forces.jl")
include("DOFevolvers.jl")
include("FieldUpdaters.jl")
include("LivePlottingFunctions.jl")


function periodic!(p_i, systemsizes)

    for (i, xi) in pairs(p_i.x)

        if xi<-systemsizes[i]/2
            p_i.x[i] = xi + systemsizes[i]
        elseif xi>=systemsizes[i]/2
            p_i.x[i] = xi - systemsizes[i]
        end
    end
    return p_i
end



function minimal_image_closest_bin_center!(field_indices,x, bin_centers,system_sizes,system_Periodic)

    
    for (i, xi) in pairs(x)
        d = @MVector zeros(length(x))
        d= abs.(bin_centers[i] .- x[i])
        # if system_Periodic
        #     d.=d .% system_sizes[i]
        # end
        field_indices[i] = findmin(d)[2]
        
    end
    return field_indices

end

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
struct System{T1, T2, T3, T4, T5, T6, T7}

    #Vector that determines the linear size of the system
    sizes::Vector{Float64}

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

struct SIM{T1, T2, T3, T4, T5}
    particle_states::T1
    field_states::T2
    tsax::T3
    dt::T4
    t_stop::T5
    system::System
end


function save_SIM(folder_path, file_name, sim)


    mkpath(folder_path)
    file_path = folder_path*file_name*".jld2"

    println(file_path)
    jldsave(file_path, true; sim=sim)
    return file_path
end

function load_SIM(file_path)
    sim_file = jldopen(file_path, "r");
    sim = sim_file["sim"]
    close(sim_file)
    return sim
end


function Euler_integrator(system, dt, t_stop,  Tsave, Tplot=nothing, fps=nothing, plot_functions=nothing,plotdim=nothing)

    
    system, cells,cell_bin_centers,stencils = construct_cell_lists!(system)

    particle_states = [copy(system.initial_particle_state)]
    field_states = [copy(system.initial_field_state)]
    tsax = [0.]


    Npair = length(system.pair_forces)

    Nfield = length(system.field_forces)

    Nfieldu = length(system.field_updaters)

    

    current_particle_state = copy(system.initial_particle_state)
    new_particle_state  = copy(system.initial_particle_state)

    current_field_state = copy(system.initial_field_state)
    new_field_state  = copy(system.initial_field_state)

    if !isnothing(plot_functions)
        if fps!=0
            cpsO = Observable(current_particle_state)
            cfsO = Observable(current_field_state)
            tO = Observable(0.)
            f, ax = setup_system_plotting(system.sizes,plot_functions, plotdim,cpsO,cfsO,tO)
        end
    end

    #Loop over time
    for (n, t) in pairs(0:dt:t_stop)

        Threads.@threads for i in eachindex(current_particle_state)
            p_i = new_particle_state[i]

            p_i, cells,new_field_state = particle_step!(i,p_i, current_particle_state,current_field_state, new_field_state,Npair, Nfield,t, dt, system,cells,cell_bin_centers,stencils)
            new_particle_state[i]=p_i
        end

        #Updating the cell list is not threadsafe, hence put it outside the threaded loop
        for i in eachindex(new_particle_state)
            p_i = new_particle_state[i]
            p_i, cells=update_cells!(p_i, cells,cell_bin_centers, stencils,system)
            cells = update_ghost_cells!(cells,system)
            new_particle_state[i]=p_i
        end
    #Now update fields
        for i in eachindex(current_field_state)

            field_i = new_field_state[i]

            field_i = field_step!(i, field_i,current_field_state,new_field_state,Nfieldu,t,dt, system)


            new_field_state[i] = field_i
                
        end

        current_particle_state = new_particle_state    
        current_field_state = new_field_state

        save_state!(particle_states, current_particle_state, n, Tsave)
        save_state!(field_states, current_field_state, n, Tsave)
        save_state!(tsax,t, n, Tsave)
        
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
    return SIM(particle_states, field_states, tsax,dt, t_stop, system)
end
function particle_step!(i,p_i, current_particle_state,current_field_state, new_field_state,Npair, Nfield,t, dt, system,cells,cell_bin_centers,stencils)
    if Npair>0
        p_i=contribute_pair_forces!(i,p_i, current_particle_state,t, dt, system,cells,stencils)
    end

    if Nfield>0
        p_i, new_field_state = contribute_field_forces!(p_i, current_field_state, new_field_state, t, dt,system)

    end

    for force in system.external_forces
        p_i=contribute_external_force!(p_i, t, dt, force)
    end

    for dofevolver in system.dofevolvers
        p_i=dofevolver(p_i, t, dt)
    end

    if system.Periodic
        p_i=periodic!(p_i, system.sizes)
    end

    return p_i,cells, new_field_state
end

function field_step!(i, field_i,current_field_state,new_field_state,Nfieldu,t,dt, system)
    if Nfieldu>0
        for field_updater in system.field_updaters
            field_i=contribute_field_update!(field_i, t, dt, field_updater)
        end
    end
    for dofevolver in system.dofevolvers
        field_i=dofevolver(field_i, t, dt)
    end
    return field_i

end
function save_state!(states, current_state, n, Tsave)

    if n%Tsave==0 && n>0
        states=push!(states,deepcopy(current_state))

    end
    return states

end

function contribute_field_forces!(p_i, current_field_state,new_field_state, t, dt,system)
    field_indices = @MVector zeros(Int,length(p_i.x))
    

    for force in system.field_forces

        for (j, field_j) in pairs(current_field_state)
            
            field_indices = minimal_image_closest_bin_center!(field_indices, p_i.x, field_j.bin_centers, system.sizes, system.Periodic)

            p_i, new_field_state[j] = contribute_field_force!(p_i, field_j, field_indices, t, dt, force)
        end
    end
    return p_i, new_field_state
end




function contribute_pair_forces!(i,p_i, current_particle_state, t, dt,system,cells,stencils)
    
    dx = @MVector zeros(Float64,length(p_i.x))
    for force in system.pair_forces

        neighbours = get_neighbours(p_i,cells,stencils)
        if !isempty(neighbours)


            for n in neighbours
                p_j = current_particle_state[n]

                if p_i.id!=p_j.id

                    dx = minimal_image_difference!(dx, p_i.x, p_j.x, system.sizes, system.Periodic)

                    dxn = norm(dx)
                    
                    if dxn<=system.rcut_pair_global
                        p_i=contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force)
                    end

                end
            end
        end
    end
    return p_i
end


function setup_system_plotting(system_sizes,plot_functions,plotdim ,cpsO,cfsO,tO,res=nothing)
    GLMakie.activate!()
    if !isnothing(res)
        f = Figure(resolution=res)
    else
        f=Figure()
    end
    title = @lift("t = $($tO)")

    if !isnothing(plotdim)
        plotdim_set = plotdim
    else
        plotdim_set = length(system_sizes)

    end

    if plotdim_set==2
        ax = Axis(f[1, 1], xlabel = "x", ylabel="y",  aspect =system_sizes[1]/system_sizes[2], title=title )
        xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax,  -system_sizes[2]/2, system_sizes[2]/2)

    elseif plotdim_set==3

        ax = Axis3(f[1, 1], xlabel = "x", ylabel="y", zlabel="z",  aspect = (1,system_sizes[2]/system_sizes[1],system_sizes[3]/system_sizes[1]), title=title)
        xlims!(ax,  -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
        zlims!(ax,  -system_sizes[3]/2, system_sizes[3]/2)

    end
    
    for plot_function in plot_functions
       ax=plot_function(ax,cpsO,cfsO)
    end
    display(f)
    return f, ax
end



function make_movie(SIM, save_path, plot_functions, fps,plotdim=nothing)

    mkpath(save_path)

    cpsO = Observable(SIM.particle_states[1])

    cfsO = Observable(SIM.field_states[1])

    tO = Observable(SIM.tsax[1])

    t_indices = range(1,length(SIM.tsax))

    fig, ax = setup_system_plotting(SIM.system.sizes,plot_functions,plotdim ,cpsO,cfsO,tO,(500,500))
    
    record(fig, save_path, t_indices; framerate=fps, compression=30) do t_index 
        cpsO[] = SIM.particle_states[t_index]
        cfsO[] = SIM.field_states[t_index]
        tO[] = SIM.tsax[t_index]
    end


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
    cells = update_ghost_cells!(cells,system)
    return system, cells, cell_bin_centers, stencils
end

function get_neighbours(p_i, cells, stencils)

    neighbours=Int64[]
    candidate=@MVector zeros(Int64, length(p_i.ci))
    for stencil in stencils

        for i in eachindex(candidate)
            candidate[i]= p_i.ci[i] + stencil[i]
        end

        if !isempty(cells[candidate...])
            for id in cells[candidate...]
                push!(neighbours,id)
            end
        end
    end
    return neighbours
end

function update_cells!(p_i, cells, cell_bin_centers, stencils,system)

    new_bin_location=@MVector zeros(Int,length(p_i.x))
    new_bin_location = minimal_image_closest_bin_center!(new_bin_location,p_i.x, cell_bin_centers,system.sizes,system.Periodic)
    if p_i.id in cells[new_bin_location...]
    else
        #remove
        filter!(e->e≠p_i.id,cells[p_i.ci...])
        #add to correct lists
        p_i.ci.= new_bin_location
        cells[new_bin_location...] = append!(cells[new_bin_location...],p_i.id)
    end
    return p_i,cells
end

function update_ghost_cells!(cells,system)
    
    if system.Periodic
        dims = length(system.sizes)
        if dims==2
            cells[1,:]=@view cells[end-1,:]
            cells[end,:]=@view cells[2,:]

            cells[:,1]=@view cells[:,end-1]
            cells[:,end]=@view cells[:,2]
        end

        if dims==3
            cells[1,:,:]=@view cells[end-1,:,:]
            cells[end,:,:]=@view cells[2,:,:]

            cells[:,1,:]=@view cells[:,end-1,:]
            cells[:,end,:]=@view cells[:,2,:]

            cells[:,:,1]=@view cells[:,:,end-1]
            cells[:,:,end]=@view cells[:,:,2]
        end


        
    end
    return cells
end
