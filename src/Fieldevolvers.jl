
function overdamped_field_evolver!(field::FuelField2d, t, dt)

    field.C.+= field.Cv*dt

    field.Cv.= field.Cf

    #reinitalize
    field.Cf.*= 0.
    return field
end
