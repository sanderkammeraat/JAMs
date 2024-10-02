
function overdamped_x_evolver!(p_i, t, dt)

    p_i.x.+= p_i.v * dt

end
function overdamped_v_evolver!(p_i, t, dt)

    p_i.v.= p_i.f/p_i.zeta

end

function overdamped_θ_evolver!(p_i, t, dt)

    p_i.θ.+= p_i.ω * dt

end 

function overdamped_ω_evolver!(p_i, t, dt)

    p_i.ω.*= 0

end 

function overdamped_f_evolver!(p_i, t, dt)

    p_i.f.*= 0

end 