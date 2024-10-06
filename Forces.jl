

using Random, Distributions
using LinearAlgebra

function minimal_image_difference(xi, xj, system_sizes, system_Periodic)

    dr = (xj - xi)
    if system_Periodic
        for (n, dx) in collect(pairs(dr))
            if dx>system_sizes[n]
                dr[n]-=system_sizes[n]
            end
            if dx<=-system_sizes[n]
                dr[n]+=system_sizes[n]
            end
        end 
    end
    return dr
end

function contribute_2d_ABP_propulsion_force!(p_i,t, dt, params, system_sizes, system_Periodic)

    p_i.f[1]+= p_i.zeta * p_i.v0 *cos(p_i.θ[1])
    p_i.f[2]+= p_i.zeta * p_i.v0 *sin(p_i.θ[1])
    #p_i.fact.+= f

    return p_i

end

function contribute_2d_ABP_angular_noise!(p_i,t, dt, params, system_sizes, system_Periodic)

    ω=sqrt(2*p_i.Dr)*rand(Normal(0, 1))

    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.ω.+= ω*sqrt(dt)/dt

    return p_i
end


function contribute_soft_disk_force!(p_i,p_j,t, dt, params, system_sizes, system_Periodic)
    
    dx = minimal_image_difference(p_i.x, p_j.x, system_sizes, system_Periodic)
    dxn = norm(dx)
    d2a = p_i.a+p_j.a
    if dxn < d2a
    f = p_i.k * (dxn-d2a) * dx/dxn
    p_i.f.+= f
     #p_i.fpas.+= f
    end
    return p_i
end



function contribute_swarm_pos_force!(p_i,p_j,t, dt, params, system_sizes, system_Periodic)

    dx = minimal_image_difference(p_i.x, p_j.x, system_sizes, system_Periodic)
    dxn = norm(dx)

    f = 1/params["N"] * (dx/dxn * (1 + params["J"]*cos.(p_j.θ-p_i.θ)[1] ) - dx/dxn^2)
    p_i.f.+= f   
    return p_i 
end

function contribute_swarm_angular_force!(p_i,p_j,t, dt, params, system_sizes, system_Periodic)

    dx = minimal_image_difference(p_i.x, p_j.x, system_sizes, system_Periodic)
    dxn = norm(dx)

    ω = 1/params["N"] * params["K"] * sin.(p_j.θ-p_i.θ)[1]/dxn
    p_i.ω.+= ω
    return p_i
end
