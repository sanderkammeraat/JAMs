
using StaticArrays

# θ for polar angle in 2d, ω for angular velocity

# p for polarity vector in 3d, q for time derivative of p

# The type field in the particle definitions must be set in increasing order, starting from 1.
abstract type Particle end

abstract type RigidBody end

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
@kwdef struct PolarParticle3d<:Particle

    id::MVector{1,Int64}
    type::MVector{1,Int64}
    m::MVector{1,Float64}   = @MVector [1.0]
    zeta::MVector{1,Float64} = @MVector [1.0]
    R::MVector{1,Float64}   = @MVector [1.0]
    
    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}

    x::MVector{3,Float64}
    xuw::MVector{3,Float64} = @MVector [0.0,0.0,0.0]
    v::MVector{3,Float64}   = @MVector [0.0,0.0,0.0]
    a::MVector{3,Float64}   = @MVector [0.0,0.0,0.0]

    f::MVector{3,Float64} = @MVector [0.0,0.0,0.0]#f->a->v->r

    p::MVector{3,Float64}
    q::MVector{3,Float64} = @MVector [0.0,0.0,0.0]
    ci::MVector{3,Int64} = @MVector [0,0,0]


end


struct PolarPolymerParticle3d<:Particle

    id::MVector{1,Int64}
    type::MVector{1,Int64}

    pol_id::MVector{1,Int64}

    id_in_pol::MVector{1, Int64}

    #Number of particles in the polymer
    pol_N::MVector{1, Int64}

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

#Parametrize on the shape extent (point defining the contour) 
struct PolarShape<:RigidBody

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
    # Updated positions of subparticles
    xe::Matrix{Float64}
    # Initial positions of subparticles relative to com
    xo::Matrix{Float64}

    # Radii of shape points for soft potential
    re::Vector{Float64}

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


#Thanks to Nina Cankilic, Chinmay Pabshettiwar and Silke Henkes
#[Reference to MSc thesis to be inserted]
struct GBEllipses<:Particle
    id::MVector{1,Int64}
    type::MVector{1,Int64}

    # Shape of ellipse
    # first element represents major axis part
    Lambda::MMatrix{2,2,Float64}
    Lambda0::MMatrix{2,2,Float64}
    sigma_0::MVector{1,Float64}

    K::MVector{1,Float64}
    mu::MVector{1,Float64}
    tau::MVector{1,Float64}

    v0::MVector{1,Float64}
    Dr::MVector{1,Float64}
    
    D::MVector{1,Float64}

    deformable::MVector{1,Bool}

    x::MVector{3,Float64}
    xuw::MVector{3,Float64}
    v::MVector{3,Float64}
    f::MVector{3,Float64} 
    stress::MMatrix{2,2,Float64}

    p::MVector{3,Float64}
    q::MVector{3,Float64}

    zeta::MVector{1,Float64} 
    fact::MVector{3,Float64}
    fpas::MVector{3,Float64}
    ci::MVector{3,Int64}
end