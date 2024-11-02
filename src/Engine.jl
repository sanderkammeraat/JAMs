
using ProgressBars
using StaticArrays

#Note, only arrays can be changed in a struct. So initializing a struct attribute as array allows to change
#Type declaration in structs is important for performance, see https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations
include("Particles.jl")
include("Forces.jl")
include("DOFevolvers.jl")
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
struct System{T1, T2, T3, T4}

    #Vector that determines the linear size of the system
    sizes::Vector{Float64}

    #Array containing particles in a specific state
    initial_state::T1

    #Array of force functions:
    external_forces::T2

    pair_forces::T3
    
    #Array of functions to evolve dof (and reinitialize forces)
    dofevolvers::T4

    #Spatially periodic boundary conditions?
    Periodic::Bool

end



function Euler_integrator(system, dt, t_stop,  Tsave, Tplot=nothing, plot_functions=nothing)

    states = [copy(system.initial_state)]

    external_forces = system.external_forces
    pair_forces = system.pair_forces

    Npair = length(pair_forces)


    current_state = copy(system.initial_state)
    new_state  = copy(system.initial_state)

    if !isnothing(plot_functions)
        if Tplot!=0
            
            f, ax = setup_plotting(system.sizes)
        end
    end

    #Loop over time
    for (n, t) in ProgressBar(pairs(0:dt:t_stop))
        #Looping over old states, so could be parallelized
        #Threads.@threads
        Np = length(current_state)
        Threads.@threads for i in 1:Np
            p_i = new_state[i]

            if Npair>0
                p_i=contribute_pair_forces!(i,p_i, current_state, pair_forces, t, dt,system.sizes, system.Periodic)
            end

            for force in external_forces
                p_i=contribute_external_force!(p_i, t, dt, force)
            end

            for dofevolver in system.dofevolvers
                p_i=dofevolver(p_i, t, dt)
            end

            if system.Periodic
                p_i=periodic!(p_i, system.sizes)
            end
            new_state[i]=p_i
        end

        current_state = new_state
        save_state!(states, current_state, n, Tsave)
        if !isnothing(plot_functions)
            if Tplot!=0
                if n%Tplot==0
                    empty!(ax)
                    

                    for plot_function in plot_functions
                        plot_function(ax, current_state)
                    end
                    ax.title="t = $(t)"
                    display(f)
                    
                end
            end
        end
    end
    return states
end

function save_state!(states, current_state, n, Tsave)

    if n%Tsave==0 && n>0
        push!(states,current_state)

    end

end


function contribute_pair_forces!(i,p_i, current_state, pair_forces, t, dt,system_sizes, system_Periodic)
    
    dx = @MVector zeros(Float64,length(p_i.x))
    for force in pair_forces
        
        for j in eachindex(current_state)

            if i!=j
                p_j = current_state[j]

                dx = minimal_image_difference!(dx, p_i.x, p_j.x, system_sizes, system_Periodic)

                dxn = norm(dx)
                
                p_i=contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force)
            end
        end
    end
    return p_i
end


function setup_plotting(system_sizes)
    GLMakie.activate!()
    f = Figure()

    dimension = length(system_sizes)
    if dimension==2
        ax = Axis(f[1, 1], xlabel = "x", ylabel="y",  aspect = 1)
        xlims!(ax, -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax,  -system_sizes[2]/2, system_sizes[2]/2)

    elseif dimension==3

        ax = Axis3(f[1, 1], xlabel = "x", ylabel="y", zlabel="z",  aspect = (1,1,1))
        xlims!(ax,  -system_sizes[1]/2, system_sizes[1]/2)
        ylims!(ax, -system_sizes[2]/2, system_sizes[2]/2)
        zlims!(ax,  -system_sizes[3]/2, system_sizes[3]/2)

    end
    display(f)

    return f, ax
end

