
using StaticArrays

# θ for polar angle in 2d, ω angular velocity

# p for polarity vector in 3d, q angular velocity

abstract type Particle end

#Defining a particle requires a class with particle properties
struct PolarParticle2d<:Particle
    id::Integer
    m::Float64 #mass, only used icw intertial integrators
    v0::Float64
    Dr::Float64

    x::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    k::Float64
    a::Float64

    zeta::Float64 #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
end

struct PolarParticle2dNtype<:Particle
    id::Integer
    type::Integer
    a::Float64
    k::Float64

    v0::Float64
    Dr::Float64

    x::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    #Neighbour memory
    n::MVector{1,Int64}
    fn::MVector{2,Float64}
    

    zeta::Float64 #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
end

struct PolarParticle3dN<:Particle
    id::Integer
    #particle radius
    a::Float64
    k::Float64
    v0::Float64
    Dr::Float64

    x::MVector{3,Float64}
    v::MVector{3,Float64}
    f::MVector{3,Float64} #f->v->r

    p::MVector{3,Float64}
    q::MVector{3,Float64}

    #Neighbour memory
    n::MVector{1,Int64}
    qn::MVector{3,Float64}
    
    

    zeta::Float64 #Friction coefficient, only used icw overdamped integrators
    fact::MVector{3,Float64}
    fpas::MVector{3,Float64}
end

struct PolarParticle3dNtype<:Particle
    id::Integer
    type::Integer
    #particle radius
    a::Float64
    k::Float64
    v0::Float64
    Dr::Float64

    x::MVector{3,Float64}
    v::MVector{3,Float64}
    f::MVector{3,Float64} #f->v->r

    p::MVector{3,Float64}
    q::MVector{3,Float64}

    #Neighbour memory
    n::MVector{1,Int64}
    qn::MVector{3,Float64}
    
    

    zeta::Float64 #Friction coefficient, only used icw overdamped integrators
    fact::MVector{3,Float64}
    fpas::MVector{3,Float64}
end


#And a way to copy this type of particle
#Here we see how the function-first approach of Julia works. We have to add to the base copy method how to act on this type (particle class)
#Base.copy(p::PolarParticle2d) = PolarParticle2d(copy(p.id), copy(p.m), copy(p.v0), copy(p.Dr), copy(p.x), copy(p.v), copy(p.f), copy(p.θ), copy(p.ω), copy(p.zeta), copy(p.fact), copy(p.fpas))

struct Swarmalator<:Particle
    id::Integer

    v0::Float64
    Dr::Float64

    x::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    a::Float64
    k::Float64


    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    ϕ::MVector{1,Float64}
    ψ::MVector{1,Float64}

    zeta::Float64 #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
end


struct VicsekParticle<:Particle
    id::Integer
    v0::Float64
    Dr::Float64

    x::MVector{2,Float64}
    v::MVector{2,Float64}
    f::MVector{2,Float64} #f->v->r

    θ::MVector{1,Float64}
    ω::MVector{1,Float64}

    #Neighbour memory
    n::MVector{1,Int64}
    ωn::MVector{1,Float64}

    zeta::Float64 #Friction coefficient, only used icw overdamped integrators
    fact::MVector{2,Float64}
    fpas::MVector{2,Float64}
end