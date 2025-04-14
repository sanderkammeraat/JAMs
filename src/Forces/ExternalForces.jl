

#Single-particle forces: External

struct external_harmonic_pinning_force<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    l::Float64
    pins::Matrix{Float64}
    
end

struct electrode_force<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    Emag::Float64
    d::Float64
    Qreset::Float64
    
end
struct self_align_with_v_force<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    β::Float64
    
end

struct self_align_with_v_unit_force<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    β::Float64
    
end

struct external_harmonic_force<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    
end

struct external_friction_force<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    γ::Float64
end

struct thermal_translational_noise<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    Dt_vector::MVector{3,Float64}
end

struct ABP_perpendicular_angular_noise<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
    perpendicular_vector::MVector{3,Float64}
end

struct ABP_2d_propulsion_force <:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
end

struct ABP_2d_angular_noise<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
end

struct ABP_3d_propulsion_force <:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
end

struct ABP_3d_angular_noise<:ExternalForce
    ontypes::Union{Int64,Vector{Int64}}
end

#Template by defining a contribute function for general type

function contribute_external_force!(p_i,t, dt, force,rngs_particles)

    if p_i.type[1] in force.ontypes
        #do force calculation
    end
    return p_i
end


function contribute_external_force!(p_i,t, dt, force::ABP_2d_propulsion_force,rngs_particles)

    if p_i.type[1] in force.ontypes
    p_i.f[1]+= p_i.zeta[1] * p_i.v0[1] *cos(p_i.θ[1])
    p_i.f[2]+= p_i.zeta[1] * p_i.v0[1] *sin(p_i.θ[1])
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::ABP_2d_angular_noise,rngs_particles)

    if p_i.type[1] in force.ontypes
    ω=sqrt(2*p_i.Dr[1])*rand(rngs_particles[p_i.id[1]],Normal(0, 1))

    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.ω.+= ω*sqrt(dt)/dt
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::thermal_translational_noise,rngs_particles)

    if p_i.type[1] in force.ontypes
        for i in eachindex(p_i.f)
            η = sqrt(2*force.Dt_vector[i])*rand(rngs_particles[p_i.id[1]], Normal(0, 1))
            #compensate for the dt from the dof evolver, can be changed if the evolver also changes
            p_i.f[i]+= η*sqrt(dt)/dt
        end
    end
    return p_i
end

function contribute_external_force!(p_i,t, dt, force::ABP_3d_propulsion_force,rngs_particles)

    if p_i.type[1] in force.ontypes
    p_i.f.+= p_i.zeta[1] * p_i.v0[1] * p_i.p
    end

    return p_i
end

function contribute_external_force!(p_i, t, dt, force::ABP_3d_angular_noise,rngs_particles)
    if p_i.type[1] in force.ontypes
    xi=sqrt(2*p_i.Dr[1])*normalize(rand(rngs_particles[p_i.id[1]],Normal(0, 1),3))

    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.q.+= cross(p_i.p, xi ) * sqrt(dt)/dt
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::ABP_perpendicular_angular_noise,rngs_particles)
    if p_i.type[1] in force.ontypes
    η =sqrt(2*p_i.Dr[1])*rand(rngs_particles[p_i.id[1]],Normal(0, 1))

    p_i.q.+= η*cross(p_i.p, force.perpendicular_vector ) * sqrt(dt)/dt
    end
    return p_i
end



function contribute_external_force!(p_i, t, dt, force::self_align_with_v_force,rngs_particles)
    if p_i.type[1] in force.ontypes
    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.q.+= force.β*cross(cross(p_i.p, p_i.v ), p_i.p)
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::self_align_with_v_unit_force,rngs_particles)
    if p_i.type[1] in force.ontypes
    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
        vnorm = norm(p_i.v)
        if vnorm!=0
            p_i.q.+= force.β*cross(cross(p_i.p,  p_i.v), p_i.p)./vnorm
        else
            p_i.q.+= force.β*cross(cross(p_i.p,  p_i.v), p_i.p)
        end
    end 
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::external_harmonic_force,rngs_particles)
    if p_i.type[1] in force.ontypes
    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.f.+= - force.k * p_i.x
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::electrode_force,rngs_particles)
    if p_i.type[1] in force.ontypes
        #compensate for the dt from the dof evolver, can be changed if the evolver also changes
        if p_i.x[3]+p_i.R[1]<=force.d/2 && p_i.x[3]-p_i.R[1]>=-force.d/2
            p_i.f[3]+= force.Emag * p_i.Q[1]

        elseif p_i.x[3]+p_i.R[1]>force.d/2
            p_i.Q.=1*force.Qreset
            p_i.v[3]*= -1
        elseif p_i.x[3]-p_i.R[1]<-force.d/2
            p_i.Q.=-1*force.Qreset
            p_i.v[3]*= -1
        end
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::external_friction_force,rngs_particles)
    if p_i.type[1] in force.ontypes
    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.f.+= - force.γ * p_i.v
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, force::external_harmonic_pinning_force,rngs_particles)
    if p_i.type[1] in force.ontypes
        dxp = @MVector zeros(length(p_i.x))
        dxp.=p_i.x - force.pins[p_i.id[1],:]
        dxpn= norm(dxp)
        if dxpn>0

            p_i.f.+= - force.k * (dxpn-force.l) * dxp/dxpn
        end
    end
    return p_i
end