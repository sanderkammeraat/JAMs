

using Random, Distributions
using LinearAlgebra
using StaticArrays


abstract type Force end


#Single-particle forces
struct field_propulsion_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    v0offset::Float64
    
end
struct field_propulsion_3d_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    v0offset::Float64
end
struct field_propulsion_distr_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    cmid::Float64
    v0max::Float64
    σ::Float64
end

struct asymmetric_field_propulsion_distr_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    cmid::Float64
    v0max::Float64
    σ::Float64
end
struct grad_field_propulsion_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    μ::Float64
end

struct self_align_with_∇C_unit_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    β::Float64
end

struct self_align_with_∇C_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    β::Float64
end

struct external_harmonic_pinning_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    l::Float64
    pins::Matrix{Float64}
    
end

struct electrode_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    Emag::Float64
    d::Float64
    Qreset::Float64
    
end
struct self_align_with_v_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    β::Float64
    
end



struct self_align_with_v_unit_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    β::Float64
    
end

struct external_harmonic_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    
end

struct external_friction_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    γ::Float64
end

struct thermal_translational_noise<:Force
    ontypes::Union{Int64,Vector{Int64}}
    Dt_vector::MVector{3,Float64}
end

struct ABP_perpendicular_angular_noise<:Force
    ontypes::Union{Int64,Vector{Int64}}
    perpendicular_vector::MVector{3,Float64}
    
end

struct ABP_2d_propulsion_force <:Force
    ontypes::Union{Int64,Vector{Int64}}
end

struct ABP_2d_angular_noise<:Force
    ontypes::Union{Int64,Vector{Int64}}
end

struct ABP_3d_propulsion_force <:Force
    ontypes::Union{Int64,Vector{Int64}}
end
struct oscillatory_3d_propulsion_force <:Force
    ontypes::Union{Int64,Vector{Int64}}
    ω::Float64
end

struct ABP_3d_angular_noise<:Force
    ontypes::Union{Int64,Vector{Int64}}
end


#Pair forces
struct coulomb_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    
end



struct soft_atre_type_force{T1, T2}
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
    ϵarray::T2
    
end


#Thanks to Julia's indexing system. karray can be a float or a 1-element vector if using only one type: type= 1
# cf. b=[2] then b[1,1] = 2 or b=2 then also b[1,1]=2. Or karray is 2d array with different stifness between different types
struct soft_disk_force{T1} <: Force
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end
struct soft_shape_disk_force{T1} <: Force
    ontypes::Union{Int64,Vector{Int64}}
    karray::T1
end

#Based on Wensink et al. PHYSICAL REVIEW E 89, 010302(R) (2014)
struct exp_shape_disk_force{T1} <: Force
    ontypes::Union{Int64,Vector{Int64}}
    u0array::T1
end



struct morse_force{T1, T2}<:Force
    ontypes::Union{Int64,Vector{Int64}}
    Dearray::T1
    aarray::T2
end


struct swarm_pos_force <: Force
    ontypes::Union{Int64,Vector{Int64}}
    N_inv::Float64
    J::Float64
    

end

struct swarm_angular_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    N_inv::Float64
    K::Float64
    

end

struct Vicsek_align_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    r::Float64
    
end

struct pairABP_force{T1}<:Force
    ontypes::Union{Int64,Vector{Int64}}
    rfact::Float64
    marray::T1
    torque::Bool
end


struct pair_polar_alignment_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    rcut::Float64
    J::Float64
end

struct chain_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    l::Float64
    
end
struct tangential_propulsion_force<:Force
    ontypes::Union{Int64,Vector{Int64}}    
end


struct spring_network_2d_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    l_network::Matrix{Float64}
    k_network::Matrix{Float64}
    
end

struct periodic_chain_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    k::Float64
    l::Float64
    periodic_id_begin::Int64
    periodic_id_end::Int64
    
end

struct fluid_dipole_2d_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    α::Float64
    l::Float64

end

struct fluid_dipole_3d_force<:Force
    ontypes::Union{Int64,Vector{Int64}}
    α::Float64
    l::Float64

end



#Let's test the power of multiple dispatch


function contribute_external_force!(p_i,t, dt, rngs_particles, system, force::ABP_2d_propulsion_force,)

    if p_i.type[1] in force.ontypes
    p_i.f[1]+= p_i.zeta[1] * p_i.v0[1] *cos(p_i.θ[1])
    p_i.f[2]+= p_i.zeta[1] * p_i.v0[1] *sin(p_i.θ[1])
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, rngs_particles, system, force::ABP_2d_angular_noise)

    if p_i.type[1] in force.ontypes
    ω=sqrt(2*p_i.Dr[1])*rand(rngs_particles[p_i.id[1]],Normal(0, 1))

    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.ω.+= ω*sqrt(dt)/dt
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt,rngs_particles,system, force::thermal_translational_noise)

    if p_i.type[1] in force.ontypes
        for i in eachindex(p_i.f)
            η = sqrt(2*force.Dt_vector[i])*rand(rngs_particles[p_i.id[1]], Normal(0, 1))
            #compensate for the dt from the dof evolver, can be changed if the evolver also changes
            p_i.f[i]+= η*sqrt(dt)/dt
        end
    end
    return p_i
end

function contribute_external_force!(p_i,t, dt,rngs_particles, system, force::ABP_3d_propulsion_force)

    if p_i.type[1] in force.ontypes
    p_i.f.+= p_i.zeta[1] * p_i.v0[1] * p_i.p
    end

    return p_i
end

function contribute_external_force!(p_i,t, dt,rngs_particles,system, force::oscillatory_3d_propulsion_force)

    if p_i.type[1] in force.ontypes
    p_i.f.+= p_i.zeta[1] * p_i.v0[1] * p_i.p * cos(force.ω*t)
    end

    return p_i
end

function contribute_external_force!(p_i, t, dt,rngs_particles,system, force::ABP_3d_angular_noise)
    if p_i.type[1] in force.ontypes
    xi=sqrt(2*p_i.Dr[1])*normalize(rand(rngs_particles[p_i.id[1]],Normal(0, 1),3))

    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.q.+= cross(p_i.p, xi ) * sqrt(dt)/dt
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt, rngs_particles, system, force::ABP_perpendicular_angular_noise)
    if p_i.type[1] in force.ontypes
    η =sqrt(2*p_i.Dr[1])*rand(rngs_particles[p_i.id[1]],Normal(0, 1))

    p_i.q.+= η*cross(p_i.p, force.perpendicular_vector ) * sqrt(dt)/dt
    end
    return p_i
end



function contribute_external_force!(p_i, t, dt,rngs_particles,system, force::self_align_with_v_force)
    if p_i.type[1] in force.ontypes
    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.q.+= force.β*cross(cross(p_i.p, p_i.v ), p_i.p)
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt,rngs_particles, system, force::self_align_with_v_unit_force)
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

function contribute_external_force!(p_i, t, dt,rngs_particles,system, force::external_harmonic_force)
    if p_i.type[1] in force.ontypes
    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.f.+= - force.k * p_i.x
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt,rngs_particles, system, force::electrode_force)
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

function contribute_external_force!(p_i, t, dt,rngs_particles, system, force::external_friction_force)
    if p_i.type[1] in force.ontypes
    #compensate for the dt from the dof evolver, can be changed if the evolver also changes
    p_i.f.+= - force.γ * p_i.v
    end
    return p_i
end

function contribute_external_force!(p_i, t, dt,rngs_particles, system, force::external_harmonic_pinning_force)
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





function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::soft_disk_force)

    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
    d2R = p_i.R[1]+p_j.R[1]
    f = @MVector zeros(length(dx))
        if dxn < d2R

            @views f.= force.karray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])] * (dxn-d2R) * dx/dxn
            p_i.f.+= f
        end
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles ,system, force::soft_shape_disk_force)

    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        f = @MVector zeros(length(dx))
        T = @MVector zeros(length(dx))


        @views for m=1:size(p_i.xe)[1]

            for n=1:size(p_j.xe)[1]
                
                d2R = p_i.re[m] + p_j.re[n]


                dxe =  @MVector zeros(length(dx))

                dxe.= minimal_image_difference!(dxe, p_i.xe[m,:],p_j.xe[n,:],system.sizes, system.Periodic)

                dxen = norm(dxe)

                if dxen<d2R
                    f.= force.karray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])] * (dxen - d2R) * dxe/dxen
                    T.=  cross( cross(p_i.xe[n,:] - p_i.x , f), p_i.p) 

                    p_i.f.+= f
                    p_i.q.+= T

                end
            end
        end
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles,system, force::exp_shape_disk_force)

    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        f = @MVector zeros(length(dx))
        T = @MVector zeros(length(dx))

        @views for m=1:size(p_i.xe)[1]

            for n=1:size(p_j.xe)[1]
                
                #Note we use here the radius not the diameter
                d2R = p_i.re[m] + p_j.re[n]


                dxe =  @MVector zeros(length(dx))

                dxe.= minimal_image_difference!(dxe, p_i.xe[m,:],p_j.xe[n,:],system.sizes, system.Periodic)

                dxen = norm(dxe)

                if dxen<d2R

                    fact = -1*exp(-dxen/d2R)*(2*d2R +dxen )/(d2R * dxen^3)

                    f.= force.u0array[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])] * fact * dxe/dxen
                    T.=  cross( cross(p_i.xe[n,:] - p_i.x , f), p_i.p) 

                    p_i.f.+= f
                    p_i.q.+= T

                end
            end
        end
    end
    return p_i

end




function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::morse_force)

    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        re = p_i.R[1]+p_j.R[1]
        f = @MVector zeros(length(dx))

        a = force.aarray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]
        De = force.Dearray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]

        @views f.= 2*De*a*( exp(-a*(dxn-re)) - exp(-2*a*(dxn-re)) ) * dx/dxn
        p_i.f.+= f

    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::chain_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        if p_j.id[1] == p_i.id[1]+1 || p_j.id[1] == p_i.id[1]-1
            p_i.f.+= force.k * (dxn-force.l) * dx/dxn
        end
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::tangential_propulsion_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        if p_j.id[1] == p_i.id[1]+1
            p_i.f.+= p_i.v0[1]*p_i.zeta[1] * dx/dxn
        end
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::spring_network_2d_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        p_i.f.+= force.k_network[p_i.id[1],p_j.id[1]] * (dxn-force.l_network[p_i.id[1],p_j.id[1]]) * dx/dxn
    end
    return p_i

end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::periodic_chain_force)
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


function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, rngs_particles, system, force::soft_atre_type_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        bij = p_i.R[1]+p_j.R[1]
        f = @MVector zeros(length(dx))
        
        ϵ = force.ϵarray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]::Float64

        r1 = (1+ϵ)*bij

        r2 = (1+2*ϵ)*bij

        if dxn <= r1
            k = force.karray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]
            f.= k * (dxn-bij) * dx/dxn
            p_i.f.+= f

        elseif  r1<dxn<r2
            k = force.karray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]
            f.=  -k * (dxn-r2) * dx/dxn

            p_i.f.+= f
        else
            0
        end
    end
    return p_i

end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::coulomb_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
    p_i.f.+= -force.k * (dx/dxn^3 * p_i.Q[1] * p_j.Q[1]) 
    end  
    return p_i
end


function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::swarm_pos_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
    @views p_i.f.+= force.N_inv * (dx/dxn * (1 + force.J*cos(p_j.ϕ[1]-p_i.ϕ[1]) ) - dx/dxn^2)  
    end 
    return p_i
end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::swarm_angular_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
    p_i.ψ.+= force.N_inv * force.K * sin(p_j.ϕ[1]-p_i.ϕ[1])/dxn
    end
    return p_i

end



function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::Vicsek_align_force)
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        if dxn < force.r
            p_i.ωn[1] += p_j.θ[1]/dt
            p_i.n[1]+=1
        end
    end
    return p_i
end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::pairABP_force)
    
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        d2a = p_i.R[1]+p_j.R[1]

        r = force.rfact*d2a::Float64

        if dxn < r

            f = @MVector zeros(length(dx))
            β = 1 - dxn/r
            f.= -β*(p_j.p *p_j.v0[1] .- p_i.p* p_i.v0[1])/2 * force.marray[get_param_ind(force.ontypes,p_i.type[1]),get_param_ind(force.ontypes,p_j.type[1])]

            p_i.f.+= f
            #add torque
            if force.torque
                p_i.q.+= cross(cross(dx/dxn .*p_i.R[1], f),p_i.p)
            end
        end
    end
    return p_i

end
function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt, rngs_particles, system, force::pair_polar_alignment_force)
    
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes
        #d2a = p_i.R[1]+p_j.R[1]

        if dxn < force.rcut

            #add torque
            p_i.q.+=force.J*cross(cross(p_i.p, p_j.p),p_i.p)
        end
    end
    return p_i

end





function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::fluid_dipole_2d_force)
    
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes

        p_i.f[1]+= -force.α * cos(p_j.θ[1]) * exp(-dxn/force.l)
        p_i.f[2]+= -force.α * sin(p_j.θ[1]) * exp(-dxn/force.l)

    end
    return p_i


end

function contribute_pair_force!(p_i, p_j, dx, dxn, t, dt,rngs_particles, system, force::fluid_dipole_3d_force)
    
    if p_i.type[1] in force.ontypes && p_j.type[1] in force.ontypes

        p_i.f.+= -force.α .* p_j.p .*1/dxn

    end
    return p_i


end

function contribute_field_force!(p_i,field_j,field_indices, t, dt,rngs_particles, system, force::field_propulsion_force)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes
        x_index = field_indices[1]
        y_index = field_indices[2]
        #print(x_index)

        p_i.f[1]+= p_i.zeta[1] * (field_j.C[x_index, y_index]+force.v0offset) *cos(p_i.θ[1])
        p_i.f[2]+= p_i.zeta[1] * (field_j.C[x_index, y_index]+force.v0offset) *sin(p_i.θ[1])

        if field_j.C[x_index, y_index]>0
            field_j.Cf[x_index, y_index]+=-force.consumption
        end
    end

    return p_i, field_j

end

function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles,system, force::field_propulsion_3d_force)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes
        x_index = field_indices[1]
        y_index = field_indices[2]
        #print(x_index)

        p_i.f[1]+= p_i.zeta[1] * (field_j.C[x_index, y_index]+force.v0offset) *p_i.p[1]
        p_i.f[2]+= p_i.zeta[1] * (field_j.C[x_index, y_index]+force.v0offset) *p_i.p[2]

        if field_j.C[x_index, y_index]>0
            field_j.Cf[x_index, y_index]+=-force.consumption
        end
    end

    return p_i, field_j

end

function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles,system, force::field_propulsion_distr_force)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes
        x_index = field_indices[1]
        y_index = field_indices[2]
        #print(x_index)



        if field_j.C[x_index, y_index]>0
            v0fact = force.v0max * force.σ/(force.σ+(log10(field_j.C[x_index, y_index]))^2)
            p_i.f[1]+= p_i.zeta[1] * (v0fact) *p_i.p[1]
            p_i.f[2]+= p_i.zeta[1] * (v0fact) *p_i.p[2]
            field_j.Cf[x_index, y_index]+=-force.consumption
        else
            field_j.Cf[x_index, y_index]=0
        end
    end

    return p_i, field_j

end

function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles,system, force::asymmetric_field_propulsion_distr_force)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes
        x_index = field_indices[1]
        y_index = field_indices[2]
        #print(x_index)



        if field_j.C[x_index, y_index]>0
            v0fact = force.v0max * force.σ/(force.σ+(log10(field_j.C[x_index, y_index]))^2)
            p_i.f[1]+= p_i.zeta[1] * (v0fact) *p_i.p[1]
            p_i.f[2]+= p_i.zeta[1] * (v0fact) *p_i.p[2]

            

            field_j.Cf[x_index-round(Int64,p_i.p[1]), y_index-round(Int64,p_i.p[2])]+=-force.consumption
        else
            field_j.Cf[x_index, y_index]=0
        end
    end

    return p_i, field_j

end

function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles, system, force::self_align_with_∇C_force)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes


        x_index = field_indices[1]
        y_index = field_indices[2]

        ∇C = grad(field_j, x_index, y_index)

        p_i.q.+= force.β*cross(cross(p_i.p, ∇C),p_i.p)
    end

    return p_i, field_j
end


function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles, system, force::self_align_with_∇C_unit_force)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes


        x_index = field_indices[1]
        y_index = field_indices[2]

        ∇C = grad(field_j, x_index, y_index)

        vnorm = norm(∇C)
        if vnorm!=0
            p_i.q.+= force.β*cross(cross(p_i.p, ∇C),p_i.p)./vnorm
        else
            p_i.q.+= force.β*cross(cross(p_i.p,  ∇C), p_i.p)
        end
    end

    return p_i, field_j
end

function contribute_field_force!(p_i,field_j,field_indices, t, dt, rngs_particles, system, force::grad_field_propulsion_force)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes


        x_index = field_indices[1]
        y_index = field_indices[2]

        ∇C = grad(field_j, x_index, y_index)

        if field_j.C[x_index, y_index]>0
            v0fact = force.μ * norm(∇C)
            p_i.f[1]+= p_i.zeta[1] * (v0fact) *p_i.p[1]
            p_i.f[2]+= p_i.zeta[1] * (v0fact) *p_i.p[2]
            field_j.Cf[x_index, y_index]+=-force.consumption
        else
            field_j.Cf[x_index, y_index]=0
        end
    end
    return p_i, field_j
end

function grad(field_j, x_index, y_index)


    grad_x = (field_j.C[x_index+1, y_index] - field_j.C[x_index-1, y_index])/2

    grad_y = (field_j.C[x_index, y_index+1] - field_j.C[x_index, y_index-1])/2

    ∇C = @MVector[grad_x, grad_y,0.]
    return ∇C
end

function get_param_ind(force_types, particle_type)
    return findfirst(isequal(particle_type),force_types)
end