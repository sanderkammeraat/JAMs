
abstract type FieldUpdater end


struct PeriodicDiffusion<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    D::Float64

end

struct EdgeSet<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float64

end

struct AvgSetwoGhost<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float64

end

struct GhostSet<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
end

function contribute_field_update!(field_i, t, dt, field_updater::PeriodicDiffusion, rngs_fields)
 
    if field_i.type in field_updater.ontypes
    Cp0 = circshift(field_i.C, (1,0))
    Cm0 = circshift(field_i.C, (-1,0))
    C0p = circshift(field_i.C, (0,1))
    C0m = circshift(field_i.C, (0,-1))

    field_i.Cf.+= field_updater.D .* (Cp0 + Cm0 +C0p +C0m   - 4 * field_i.C )
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

function contribute_field_update!(field_i, t, dt, field_updater::AvgSetwoGhost, rngs_fields)

    if field_i.type in field_updater.ontypes

        avg = mean(field_i.C[2:end-1,2:end-1])
        field_i.C[2:end-1,2:end-1].+= field_updater.Cset - avg
    end
    return field_i
end


function contribute_field_update!(field_i, t, dt, field_updater::GhostSet, rngs_fields)

    if field_i.type in field_updater.ontypes

    field_i.C[1,2:end-1].= field_i.C[end-1,2:end-1]

    field_i.C[end,2:end-1].= field_i.C[2,2:end-1]

    field_i.C[2:end-1,1].= field_i.C[2:end-1,end-1]

    field_i.C[2:end-1,end].= field_i.C[2:end-1,2]

    field_i.C[1,1]= field_i.C[end-1,end-1]

    field_i.C[end,end]= field_i.C[2,2]

    field_i.C[end,1]= field_i.C[2,end-1]

    field_i.C[1,end]= field_i.C[end-1,2]


    field_i.Cf[1,2:end-1].= field_i.Cf[end-1,2:end-1]

    field_i.Cf[end,2:end-1].= field_i.Cf[2,2:end-1]

    field_i.Cf[2:end-1,1].= field_i.Cf[2:end-1,end-1]

    field_i.Cf[2:end-1,end].= field_i.Cf[2:end-1,2]

    field_i.Cf[1,1]= field_i.Cf[end-1,end-1]

    field_i.Cf[end,end]= field_i.Cf[2,2]

    field_i.Cf[end,1]= field_i.Cf[2,end-1]

    field_i.Cf[1,end]= field_i.Cf[end-1,2]



    field_i.Cv[1,2:end-1].= field_i.Cv[end-1,2:end-1]

    field_i.Cv[end,2:end-1].= field_i.Cv[2,2:end-1]

    field_i.Cv[2:end-1,1].= field_i.Cv[2:end-1,end-1]

    field_i.Cv[2:end-1,end].= field_i.Cv[2:end-1,2]

    field_i.Cv[1,1]= field_i.Cv[end-1,end-1]

    field_i.Cv[end,end]= field_i.Cv[2,2]

    field_i.Cv[end,1]= field_i.Cv[2,end-1]

    field_i.Cv[1,end]= field_i.Cv[end-1,2]





    end
    return field_i
end