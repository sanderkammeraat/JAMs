
abstract type FieldUpdater end


struct PeriodicDiffusion<:FieldUpdater
    D::Float64

end

function contribute_field_update!(field_i, t, dt, field_updater::PeriodicDiffusion)

    Cp0 = circshift(field_i.C, (1,0))
    Cm0 = circshift(field_i.C, (-1,0))
    C0p = circshift(field_i.C, (0,1))
    C0m = circshift(field_i.C, (0,-1))

    field_i.Cf.+= field_updater.D .* (Cp0 + Cm0 +C0p +C0m   - 4 * field_i.C )

    return field_i
end