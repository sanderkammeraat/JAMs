
function overdamped_evolver!(p_i::PolarParticle2d, t, dt)

    #Evolve
    p_i.x.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta
    p_i.θ.+= p_i.ω * dt

    #reinitalize
    p_i.ω.*= 0.
    p_i.f.*= 0.

    return p_i

end

function overdamped_evolver!(p_i::Swarmalator, t, dt)

    #Evolve
    p_i.x.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta
    p_i.θ.+= p_i.ω * dt

    #reinitalize
    p_i.ω.*= 0.
    p_i.f.*= 0.

    return p_i

end

function overdamped_evolver!(p_i::VicsekParticle, t, dt)

    #Evolve
    #Process neighbour memory:

    if p_i.n[1]>0
        p_i.ω[1]+= p_i.ωn[1]/p_i.n[1]
    end

    #
    p_i.x.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta
    p_i.θ.= p_i.ω * dt

    #reinitalize
    p_i.ω.*= 0.
    p_i.f.*= 0.

    p_i.ωn[1]= 0.
    p_i.n[1]= 0

    return p_i

end
