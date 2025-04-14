#Multi-particle forces: Pair
struct coulomb_force<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    
end

struct soft_atre_type_force{T1, T2}<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
    ϵarray::T2
    
end

#Thanks to Julia's indexing system. karray can be a float or a 1-element vector if using only one type: type= 1
# cf. b=[2] then b[1,1] = 2 or b=2 then also b[1,1]=2. Or karray is 2d array with different stifness between different types
struct soft_disk_force{T1} <: PairForce
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end

struct chain_force<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    l::Float64
    
end

struct spring_network_2d_force<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    l::Float64
    k_network::Matrix{Float64}
    
end

struct periodic_chain_force<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    l::Float64
    periodic_id_begin::Int64
    periodic_id_end::Int64
    
end

struct fluid_dipole_2d_force<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    α::Float64
    l::Float64

end

struct fluid_dipole_3d_force<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    α::Float64
    l::Float64

end

struct pairABP_force<:PairForce
    ontypes::Union{Int64,Vector{Int64}}
    rfact::Float64    
end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force,rngs_particles)

    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
    #Do calculation on p_i
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::soft_disk_force,rngs_particles)

    if p_i.type[1]in force.ontypes && p_j.type[1] in force.ontypes
    d2R = p_i.R[1]+p_j.R[1]
    f = @MVector zeros(length(dx))
        if dxn < d2R

            @views f.= force.karray[p_i.type[1],p_j.type[1]] * (dxn-d2R) * dx/dxn
            p_i.f.+= f
        end
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::chain_force,rngs_particles)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        if p_j.id[1] == p_i.id[1]+1 || p_j.id[1] == p_i.id[1]-1
            p_i.f.+= force.k * (dxn-force.l) * dx/dxn
        end
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::spring_network_2d_force,rngs_particles)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        p_i.f.+= force.k_network[p_i.id[1],p_j.id[1]] * (dxn-force.l) * dx/dxn
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::periodic_chain_force,rngs_particles)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        if p_j.id[1] == p_i.id[1]+1 || p_j.id[1] == p_i.id[1]-1
            p_i.f.+= force.k * (dxn-force.l) * dx/dxn
        end

        if p_i.id[1] ==force.periodic_id_begin  && p_j.id[1] ==force.periodic_id_end
            p_i.f.+= force.k * (dxn-force.l) * dx/dxn
        end

        if p_i.id[1] ==force.periodic_id_end  && p_j.id[1] ==force.periodic_id_begin
            p_i.f.+= force.k * (dxn-force.l) * dx/dxn
        end
    end
    return p_i

end


function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::soft_atre_type_force,rngs_particles)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        bij = p_i.R[1]+p_j.R[1]
        f = @MVector zeros(length(dx))
        
        ϵ = force.ϵarray[p_i.type[1],p_j.type[1]]::Float64

        r1 = (1+ϵ)*bij

        r2 = (1+2*ϵ)*bij

        if dxn <= r1
            k = force.karray[p_i.type[1],p_j.type[1]]

            f.= k * (dxn-bij) * dx/dxn
            p_i.f.+= f

        elseif  r1<dxn<r2
            k = force.karray[p_i.type[1],p_j.type[1]]
            f.=  k * (dxn-r2) * dx/dxn

            p_i.f.+= f
        else
            0
        end
    end
    return p_i

end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::coulomb_force,rngs_particles)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
    p_i.f.+= -force.k * (dx/dxn^3 * p_i.Q[1] * p_j.Q[1]) 
    end  
    return p_i
end



function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::fluid_dipole_2d_force,rngs_particles)
    
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes

        p_i.f[1]+= -force.α * cos(p_j.θ[1]) * exp(-dxn/force.l)
        p_i.f[2]+= -force.α * sin(p_j.θ[1]) * exp(-dxn/force.l)

    end
    return p_i


end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, force::fluid_dipole_3d_force,rngs_particles)
    
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes

        p_i.f.+= -force.α .* p_j.p .*1/dxn

    end
    return p_i


end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,  force::pairABP_force,rngs_particles)
    
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        d2a = p_i.R[1]+p_j.R[1]

        r = force.rfact*d2a

        if dxn < r
            β = 1 - dxn/r
            p_i.f. += -β*(p_j.p .*p_j.v0[1] - p_i.p .*p_i.v0[1])/2
        end
    end
    return p_i
end