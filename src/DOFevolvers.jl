function inertial_evolver!(p_i::Hexbug, t, dt)

    #evolve
    p_i.a .= p_i.f/p_i.m

    p_i.x .+= p_i.v * dt
    p_i.v .+= p_i.a * dt

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end
function inertial_evolver!(p_i::ChargedParticle3d, t, dt)

    #evolve
    p_i.a .= p_i.f/p_i.m

    p_i.x .+= p_i.v * dt
    p_i.v .+= p_i.a * dt

    #reinitialize
    p_i.f.*= 0.

    return p_i
end
function overdamped_evolver!(p_i::Hexbug, t, dt)

    #evolve
    p_i.x .+= p_i.v * dt
    p_i.v .= p_i.f/p_i.zeta

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end
function inertial_evolver!(p_i::NewPolarParticle3d, t, dt)

    #evolve
    p_i.a .= p_i.f/p_i.m

    p_i.x .+= p_i.v * dt
    p_i.v .+= p_i.a * dt

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end

function overdamped_evolver!(p_i::NewPolarParticle3d, t, dt)

    #evolve
    p_i.x .+= p_i.v * dt
    p_i.v .= p_i.f/p_i.zeta

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end



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

function overdamped_evolver!(p_i::NewPolarParticle3d, t, dt)

    #evolve
    p_i.x .+= p_i.v * dt
    p_i.v .= p_i.f/p_i.zeta

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end



function overdamped_evolver!(p_i::PolarParticle2dNtype, t, dt)

    if p_i.n[1]>0
        p_i.f.+= p_i.fn/p_i.n[1]
    end

    #Evolve
    p_i.x.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta
    p_i.θ.+= p_i.ω * dt

    #reinitalize
    p_i.ω.*= 0.
    p_i.f.*= 0.
    p_i.fn.*=0.
    p_i.n[1]= 0

    return p_i

end

function overdamped_evolver!(p_i::Swarmalator, t, dt)

    #Evolve
    p_i.x.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta
    p_i.θ.+= p_i.ω * dt

    p_i.ϕ.+= p_i.ψ * dt

    #reinitalize
    p_i.ω.*= 0.
    p_i.ψ.*= 0.
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


function overdamped_evolver!(p_i::PolarParticle3dN, t, dt)

    #Evolve
    #Process neighbour memory:

    if p_i.n[1]>0
        p_i.q.+= p_i.qn/p_i.n[1]
    end

    #
    p_i.x.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta
    p_i.p.+= p_i.q * dt
    p_i.p.=normalize(p_i.p)

    #reinitalize
    p_i.q.*= 0.
    p_i.f.*= 0.

    p_i.qn.*= 0.
    p_i.n[1]= 0

    return p_i

end


function overdamped_evolver!(field::FuelField2d, t, dt)

    field.C.+= field.Cv*dt

    field.Cv.= field.Cf

    #reinitalize
    field.Cf.*= 0.
    return field
end