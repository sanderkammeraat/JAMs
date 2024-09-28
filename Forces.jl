using Random, Distributions
using LinearAlgebra

function contribute_2d_ABP_propulsion_force!(p_i,t, dt)
    f = p_i.zeta * p_i.v0 * [cos(p_i.θ[1]), sin(p_i.θ[1])]
    p_i.f.+= f
    #p_i.fact.+= f
end

function contribute_2d_ABP_angular_noise!(p_i,t, dt)

    ω=sqrt(2*p_i.Dr)*rand(Normal(0, 1))

    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.ω.= ω*sqrt(dt)/dt
end

function contribute_soft_disk_force!(p_i,p_j,t, dt, system_sizes)
    
    dr = (p_j.x - p_i.x) .% system_sizes
    drn = norm(dr)
    d2a = p_i.a+p_j.a
    if drn < d2a
    f = p_i.k * (drn-d2a) * dr/drn
    p_i.f.+= f
     #p_i.fpas.+= f
    end
end
