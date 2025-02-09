
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
        if system_Periodic
            d.=d .% system_sizes[i]
        end
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


function Euler_integrator(system, dt, t_stop,  Tsave, Tplot=nothing, fps=nothing, plot_functions=nothing,plot_on_plane=false)

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
            f, ax = setup_system_plotting(system.sizes,plot_functions, plot_on_plane,cpsO,cfsO,tO)
        end
    end

    #Loop over time
    for (n, t) in pairs(0:dt:t_stop)

        Np = length(current_particle_state)
        Threads.@threads for i in 1:Np
            p_i = new_particle_state[i]

            p_i, new_field_state = particle_step!(i,p_i, current_particle_state,current_field_state, new_field_state,Npair, Nfield,t, dt, system)
            new_particle_state[i]=p_i
        end

    #Now update fields
        Nf = length(current_field_state)
        for i in 1:Nf

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
function particle_step!(i,p_i, current_particle_state,current_field_state, new_field_state,Npair, Nfield,t, dt, system)
    if Npair>0
        p_i=contribute_pair_forces!(i,p_i, current_particle_state,t, dt, system)
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
    return p_i, new_field_state
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




function contribute_pair_forces!(i,p_i, current_particle_state, t, dt,system)
    
    dx = @MVector zeros(Float64,length(p_i.x))
    for force in system.pair_forces
        
        for j in eachindex(current_particle_state)

            if i!=j
                p_j = current_particle_state[j]

                dx = minimal_image_difference!(dx, p_i.x, p_j.x, system.sizes, system.Periodic)

                dxn = norm(dx)
                
                if dxn<=system.rcut_pair_global
                    p_i=contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force)
                end

            end
        end
    end
    return p_i
end


function setup_system_plotting(system_sizes,plot_functions, plot_on_plane,cpsO,cfsO,tO)
    GLMakie.activate!()
    f = Figure()
    dimension = length(system_sizes)
    title = @lift("t = $($tO)")
    if dimension==2 || plot_on_plane
        ax = Axis(f[1, 1], xlabel = "x", ylabel="y",  aspect =system_sizes[1]/system_sizes[2], title=title )
        xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax,  -system_sizes[2]/2, system_sizes[2]/2)

    elseif dimension==3

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
