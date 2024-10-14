

using Random, Distributions
using LinearAlgebra
using StaticArrays


abstract type Force end


struct ABP_2d_propulsion_force <:Force
end

struct ABP_2d_angular_noise<:Force
end


struct soft_disk_force <: Force
end

struct swarm_pos_force <: Force

    N_inv::Float64
    J::Float64

end

struct swarm_angular_force<:Force

    N_inv::Float64
    K::Float64

end




#Let's test the power of multiple dispatch


function contribute_external_force!(p_i,t, dt, force::ABP_2d_propulsion_force)

    p_i.f[1]+= p_i.zeta * p_i.v0 *cos(p_i.θ[1])
    p_i.f[2]+= p_i.zeta * p_i.v0 *sin(p_i.θ[1])

    return p_i
end

function contribute_external_force!(p_i, t, dt, force::ABP_2d_angular_noise)

    ω=sqrt(2*p_i.Dr)*rand(Normal(0, 1))

    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.ω.+= ω*sqrt(dt)/dt

    return p_i
end



function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::soft_disk_force)

    d2a = p_i.a+p_j.a
    f = @MVector zeros(length(dx))
    if dxn < d2a

        f.= p_i.k * (dxn-d2a) * dx/dxn
        p_i.f.+= f
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::swarm_pos_force)

    p_i.f.+= force.N_inv * (dx/dxn * (1 + force.J*cos(p_j.θ[1]-p_i.θ[1]) ) - dx/dxn^2)   
    return p_i
end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::swarm_angular_force)

    p_i.ω.+= force.N_inv * force.K * sin(p_j.θ[1]-p_i.θ[1])/dxn
    return p_i

end
