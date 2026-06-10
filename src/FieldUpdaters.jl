
abstract type FieldUpdater end

#Set: directly set
# Otherwise, update the time derivative of the field

struct PeriodicDiffusion<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    D::Float64

end

struct Diffusion<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    D::Float64

end
struct GPUDiffusion<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    D::Float32

end
struct CPU_to_GPU<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
end

struct EdgeSet<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float64

end

struct IndSet{T1, T2}<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    inds::T1
    Cset::T2

end

struct AvgSetwoGhost<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float64

end

struct Relax<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float64
    k::Float64

end

struct GPURelax<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float32
    k::Float32

end


struct GhostSet<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
end

# function contribute_field_update!(field_i, t, dt, field_updater::PeriodicDiffusion, rngs_fields)
 
#     if field_i.type in field_updater.ontypes
#     Cp0 = circshift(field_i.C, (1,0))
#     Cm0 = circshift(field_i.C, (-1,0))
#     C0p = circshift(field_i.C, (0,1))
#     C0m = circshift(field_i.C, (0,-1))

#     field_i.Cf.+= field_updater.D .* (Cp0 .+ Cm0  .+ C0p .+ C0m   .- 4 .* field_i.C )
#     end
#     return field_i
# end

function contribute_field_update!(field_i, t, dt, field_updater::Diffusion, rngs_fields)
 
    if field_i.type in field_updater.ontypes

        Dfactor = field_updater.D/field_i.lbin^2

        C = field_i.C
        Cf = field_i.Cf


        Threads.@threads for j in 2:size(field_i.C,2)-1

            @inbounds for i in 2:size(field_i.C,1)-1

                Cf[i,j] +=  Dfactor * (C[i+1,j] + C[i-1,j]  + C[i,j+1] + C[i,j-1]  - 4 * C[i,j] )
            end
        end
    end
    return field_i
end



function contribute_field_update!(field_i, t, dt, field_updater::EdgeSet, rngs_fields)

    if field_i.type in field_updater.ontypes
    field_i.C[1,:].= field_updater.Cset
    field_i.Cf[1,:].= 0

    field_i.C[:,1].= field_updater.Cset
    field_i.Cf[:,1].= 0

    field_i.C[end,:].= field_updater.Cset
    field_i.Cf[end,:].= 0

    field_i.C[:,end].= field_updater.Cset
    field_i.Cf[:,end].= 0
    end
    return field_i
end


function contribute_field_update!(field_i, t, dt, field_updater::IndSet, rngs_fields)

    if field_i.type in field_updater.ontypes
    field_i.C[field_updater.inds...]= field_updater.Cset

    end
    return field_i
end


function contribute_field_update!(field_i, t, dt, field_updater::AvgSetwoGhost, rngs_fields)

    if field_i.type in field_updater.ontypes

        avg = mean(field_i.C[2:end-1,2:end-1])
        field_i.C[2:end-1,2:end-1].+= field_updater.Cset - avg
    end
    return field_i
end

function contribute_field_update!(field_i, t, dt, field_updater::Relax, rngs_fields)

    if field_i.type in field_updater.ontypes

        C = field_i.C
        Cf = field_i.Cf

            Threads.@threads for j in 2:size(field_i.C,2)-1

                @inbounds for i in 2:size(field_i.C,1)-1

                    Cf[i,j] += field_updater.k * (field_updater.Cset - C[i,j])

                end

            end
    end
    return field_i
end
@kernel function relax_kernel!(Cf, C,k, Cset)

    i, j = @index(Global, NTuple)

    if i > 1 && i < size(C, 1) && j > 1 && j < size(C, 2)
        @inbounds Cf[i, j] += k * (Cset - C[i, j]
        )
    end
end


function contribute_field_update!(field_i, t, dt, field_updater::GPURelax, rngs_fields)

    C_GPU = field_i.C_GPU
    Cf_GPU = field_i.Cf_GPU
    

    backend = KernelAbstractions.get_backend(C_GPU)

    if t==0
        display(backend)
    end

    kernel! = relax_kernel!(backend)

    kernel!(Cf_GPU, C_GPU, field_updater.k,field_updater.Cset, ndrange=size(C_GPU))
    
    KernelAbstractions.synchronize(backend)
    return field_i
end



@kernel function diffusion_kernel!(Cf, C,Dfactor)

    i, j = @index(Global, NTuple)

    if i > 1 && i < size(C, 1) && j > 1 && j < size(C, 2)
        @inbounds Cf[i, j] += Dfactor * (
            C[i+1, j] + C[i-1, j] + 
            C[i, j+1] + C[i, j-1] - 
            4 * C[i, j]
        )
    end
end


function contribute_field_update!(field_i, t, dt, field_updater::CPU_to_GPU, rngs_fields)

    if field_i.type in field_updater.ontypes
        copyto!(field_i.Cf_GPU, field_i.Cf)
    end
    return field_i
end

function contribute_field_update!(field_i, t, dt, field_updater::GPUDiffusion, rngs_fields)

    if field_i.type in field_updater.ontypes

        Dfactor = field_updater.D/field_i.lbin^2
        
        C_GPU = field_i.C_GPU
        Cf_GPU = field_i.Cf_GPU
        
        backend = KernelAbstractions.get_backend(C_GPU)

        if t==0
            display(backend)
        end

        kernel! = diffusion_kernel!(backend)


        kernel!(Cf_GPU, C_GPU, Dfactor, ndrange=size(C_GPU))
        

        KernelAbstractions.synchronize(backend)

    end
    return field_i
end

#Need a GPU version of this
function contribute_field_update!(field_i, t, dt, field_updater::GhostSet, rngs_fields)

    if field_i.type in field_updater.ontypes

        set_ghost_values!(field_i.C)

        set_ghost_values!(field_i.Cv)

        set_ghost_values!(field_i.Cf)

    end
    return field_i
end

function set_ghost_values!(C)

    #corner points and edges are correctly set, see ghost cell update comments in engine
    C[1,:] .= @view C[end-1,:]
    C[end,:] .= @view C[2,:]

    C[:,1] .= @view C[:,end-1]
    C[:,end] .= @view C[:,2]


    return C
end