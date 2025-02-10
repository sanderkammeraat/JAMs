
abstract type FieldUpdater end


struct PeriodicDiffusion<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    D::Float64

end

struct EdgeSet<:FieldUpdater
    ontypes::Union{Int64,Vector{Int64}}
    Cset::Float64

end

function contribute_field_update!(field_i, t, dt, field_updater::PeriodicDiffusion)
 
    if field_i.type in field_updater.ontypes
    Cp0 = circshift(field_i.C, (1,0))
    Cm0 = circshift(field_i.C, (-1,0))
    C0p = circshift(field_i.C, (0,1))
    C0m = circshift(field_i.C, (0,-1))

    field_i.Cf.+= field_updater.D .* (Cp0 + Cm0 +C0p +C0m   - 4 * field_i.C )
    end
    return field_i
end



function contribute_field_update!(field_i, t, dt, field_updater::EdgeSet)

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