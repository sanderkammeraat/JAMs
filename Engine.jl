using ProgressBars
using Plots
#Set default plotting backend
gr()
#Note, only arrays can be changed in a struct. So initializing a struct attribute as array allows to change
#Type declaration in structs is important for performance, see https://docs.julialang.org/en/v1/manual/performance-tips/#Type-declarations

function periodic!(p_i, systemsizes)

    for (i, xi) in pairs(p_i.x)

        if xi<0
            p_i.x[i] = xi + systemsizes[i]
        elseif xi>systemsizes[i]
            p_i.x[i] = xi - systemsizes[i]
        end
    end
    return p_i
end

struct Force

    #Name of force in the system
    #Maybe not needed?
    name::String

    #Is it a "pair", "graph" or "external" force?
    kind::String

    params::Dict{String,Float64}

    #If pair: (p_i, p_j, t) elif external (p_i, t) elif graph (p_i, graph, t)
    #Include t in argument even if t is not used for the calculation
    contribute::Function

end
struct System

    #Vector that determines the linear size of the system
    sizes::Vector{Float64}

    #Array containing particles in a specific state
    initial_state::AbstractVector

    #Array of force functions:
    forces::Vector{Force}

    #Array of functions to evolve dof (and reinitialize forces)
    dofevolvers::Vector{Function}

    #Spatially periodic boundary conditions?
    Periodic::Bool

end


function get_forces_of_kind(system, kind)
    #return an array with the forces if they are of the right kind
    return [force for force in system.forces if force.kind==kind]

end

function Euler_integrator(system, dt, t_stop,  Tsave, Tplot=nothing, plot_state=nothing)

    states = [copy(system.initial_state)]

    pair_forces = get_forces_of_kind(system, "pair")

    external_forces = get_forces_of_kind(system, "external")

    Npair = length(pair_forces)

    
    current_state = copy(system.initial_state)
    new_state  = copy(system.initial_state)

    if !isnothing(plot_state)
        if Tplot!=0
            p = plot(aspect_ratio=:equal)
            xlims!(0, system.sizes[1])
            ylims!(0, system.sizes[2])

            xlabel!("x")
            ylabel!("y")
            display(p)
        end
    end

    #Loop over time
    for (n, t) in ProgressBar(pairs(0:dt:t_stop))
        #Looping over old states, so could be parallelized
        #Threads.@threads
        Threads.@threads for i in eachindex(current_state)
            p_i = new_state[i]

            if Npair>0
    
                for j in eachindex(current_state)

                    if i!=j

                        for force in pair_forces
                            force.contribute(p_i, current_state[j], t, dt, force.params, system.sizes, system.Periodic)
                        end
                    end
                end
            end
            for force in external_forces
                force.contribute(p_i, t, dt, force.params, system.sizes, system.Periodic)
            end

            for dofevolver in system.dofevolvers
                dofevolver(p_i, t, dt)
            end

            if system.Periodic
                periodic!(p_i, system.sizes)
            end
            new_state[i]=p_i
        end
        current_state = new_state
        save_state!(states, current_state, n, Tsave)
        if !isnothing(plot_state)
            if Tplot!=0
                plot_state(p, current_state, n, Tplot)
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
