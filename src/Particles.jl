
using StaticArrays

# θ for polar angle in 2d, ω for angular velocity

# p for polarity vector in 3d, q for time derivative of p

# The type field in the particle definitions must be set in increasing order, starting from 1.
abstract type Particle end

struct Hexbug<:Particle

    id::MVector{1,Int64}
    type::MVector{1,Int64}
    m::MVector{1,Float64}
    zeta::MVector{1,Float64}
    
    k::MVector{1,Float64}
    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{3,Float64}
    xuw::MVector{3,Float64}
    v::MVector{3,Float64}
    a::MVector{3,Float64}

    f::MVector{3,Float64} #f->a->v->r

    p::MVector{3,Float64}
    q::MVector{3,Float64}
    ci::MVector{3,Int64}


end
struct ChargedParticle3d<:Particle
    id::MVector{1,Int64}
    type::MVector{1,Int64}
    m::MVector{1,Float64}
    zeta::MVector{1,Float64}
    R::MVector{1,Float64}
    Q::MVector{1,Float64}

    x::MVector{3,Float64}
    xuw::MVector{3,Float64}
    v::MVector{3,Float64}
    a::MVector{3,Float64}
    f::MVector{3,Float64} #f->a->v->r
    ci::MVector{3,Int64}

end
struct PolarParticle3d<:Particle

    id::MVector{1,Int64}
    type::MVector{1,Int64}
    m::MVector{1,Float64}
    zeta::MVector{1,Float64}
    R::MVector{1,Float64}
    
    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{3,Float64}
    xuw::MVector{3,Float64}
    v::MVector{3,Float64}
    a::MVector{3,Float64}

    f::MVector{3,Float64} #f->a->v->r

    p::MVector{3,Float64}
    q::MVector{3,Float64}
    ci::MVector{3,Int64}


end

struct ConfinedPolarParticle3d<:Particle

    id::MVector{1,Int64}
    type::MVector{1,Int64}
    m::MVector{1,Float64}
    zeta::MVector{1,Float64}
    R::MVector{1,Float64}
    
    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{3,Float64}
    xuw::MVector{3,Float64}
    v::MVector{3,Float64}
    a::MVector{3,Float64}

    f::MVector{3,Float64} #f->a->v->r

    p::MVector{3,Float64}
    q::MVector{3,Float64}
    ci::MVector{3,Int64}


end


#Defining a particle requires a class with particle properties
struct PolarParticle2d<:Particle
    id::MVector{1,Int64}
    type::MVector{1,Int64}
    m::MVector{1,Float64} #mass, only used icw intertial integrators
    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{2,Float64}
    xuw::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    R::MVector{1,Float64}

    zeta::MVector{1,Float64} #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
    ci::MVector{2,Int64}
end


struct PolarParticle2dN<:Particle
    id::MVector{1,Int64}
    type::MVector{1,Int64}
    R::MVector{1,Float64}
    k::MVector{1,Float64}

    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{2,Float64}
    xuw::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    #Neighbour memory
    n::MVector{1,Int64}
    fn::MVector{2,Float64}
    

    zeta::MVector{1,Float64} #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
    ci::MVector{2,Int64}
end

struct PolarParticle3dN<:Particle
    id::MVector{1,Int64}
    type::MVector{1,Int64}

    #particle radius
    R::MVector{1,Float64}
    k::MVector{1,Float64}
    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{3,Float64}
    xuw::MVector{3,Float64}
    v::MVector{3,Float64}
    f::MVector{3,Float64} #f->v->r

    p::MVector{3,Float64}
    q::MVector{3,Float64}

    #Neighbour memory
    n::MVector{1,Int64}
    qn::MVector{3,Float64}
    
    
    

    zeta::MVector{1,Float64} #Friction coefficient, only used icw overdamped integrators
    fact::MVector{3,Float64}
    fpas::MVector{3,Float64}
    ci::MVector{3,Int64}
end

struct Swarmalator<:Particle
    id::MVector{1,Int64}
    type::MVector{1,Int64}


    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{2,Float64}
    xuw::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    R::MVector{1,Float64}
    k::MVector{1,Float64}


    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    ϕ::MVector{1,Float64}
    ψ::MVector{1,Float64}

    zeta::MVector{1,Float64} #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
    ci::MVector{2,Int64}
end


struct VicsekParticle<:Particle
    id::MVector{1,Int64}
    type::MVector{1,Int64}

    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{2,Float64}
    xuw::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    #Neighbour memory
    n::MVector{1,Int64}
    ωn::MVector{1,Float64}

    zeta::MVector{1,Float64} #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
    ci::MVector{2,Int64}
end
