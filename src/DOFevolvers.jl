
#new, local_evolver and global_evolver

struct overdamped_pairdis_evolver
    η::Float64
    range::Float64 #should be smaller than the cut off radius
end

struct overdamped_xvf_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_pq_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_θω_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

#new, focus on dof but individual
function evolve_locally!(p_i, t, dt, dofevolver::overdamped_xvf_evolver)

    if p_i.type[1] in dofevolver.ontypes

        #evolve
        p_i.x .+= p_i.v * dt
        p_i.xuw.+= p_i.v * dt
        p_i.v .= p_i.f/p_i.zeta[1]

        #reinitialize
        p_i.f.*= 0.
    end
    return p_i
end

function evolve_locally!(p_i, t, dt, dofevolver::overdamped_pq_evolver)

    if p_i.type[1] in dofevolver.ontypes
        #evolve
        p_i.p .+= p_i.q * dt
        p_i.p .=normalize(p_i.p)

        #reinitialize
        p_i.q.*= 0.
    end
    return p_i
end

function evolve_locally!(p_i, t, dt, dofevolver::overdamped_θω_evolver)

    if p_i.type[1] in dofevolver.ontypes
    #evolve
    p_i.θ.+= p_i.ω * dt

    #reinitalize
    p_i.ω.*= 0.
    end

    return p_i
end

function evolve_globally!(current_particle_state, current_field_state, system, cells, stencils, dt, dofevolver::overdamped_pairdis_evolver)

    dims = length(system.sizes)
    Np = length(current_particle_state)
    id = SMatrix{dims, dims}(1I)
    id_global = SMatrix{dims, dims}(dofevolver.η*I) 

    #We are going to solve Mx=b with b the forces (excluding pair-dissipation) and M = I \zeta + M_pairdiss.
    # For that we need to construct M, first preallocate (make sparse?)

    M = spzeros(Float64,dims*Np, dims*Np)

    for i in eachindex(current_particle_state)
        p_i = current_particle_state[i]

        dx = @MVector zeros(Float64,length(p_i.x))
        neighbours= get_neighbours(p_i,cells,stencils)
        if !isnothing(neighbours)
    
            for n in neighbours
    
                if i!=n
                    p_j = current_particle_state[n]
    
                    dx = minimal_image_difference!(dx, p_i.x, p_j.x, system.sizes, system.Periodic)
    
                    dxn = norm(dx)
                    
                    @views if dxn<=dofevolver.range #In contact!
                        M[dims*(i-1)+1:dims*i, dims*(n-1)+1:dims*n].-=id_global

                        M[dims*(i-1)+1:dims*i, dims*(i-1)+1:dims*i].+= id_global
                    end
                end
                if i==n
                    M[dims*(i-1)+1:dims*i, dims*(i-1)+1:dims*i].+= id .*p_i.zeta[1]
                end
            end    
        end
    end
    b = reduce(vcat,[p_i.f for p_i in current_particle_state])
    vsol = Symmetric(M) \ b

    Threads.@threads for i in eachindex(current_particle_state)

        p_i = current_particle_state[i]

        p_i.v.= vsol[dims*(i-1)+1:dims*i]

        p_i.x .+= p_i.v * dt
        p_i.xuw.+= p_i.v * dt
        #reinitialize
        p_i.f.*= 0.

        current_particle_state[i] = p_i
    end
    return current_particle_state, current_field_state
end






function overdamped_pq_evolver!(p_i, t, dt)
    #evolve
    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)

    #reinitialize
    p_i.q.*= 0.
    return p_i
end



#old, still fully functional, focus on particle

function inertial_evolver!(p_i::Hexbug, t, dt)

    #evolve
    p_i.a .= p_i.f/p_i.m[1]

    p_i.x .+= p_i.v * dt
    p_i.xuw.+= p_i.v * dt

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
    p_i.a .= p_i.f/p_i.m[1]

    p_i.x .+= p_i.v * dt
    p_i.xuw .+= p_i.v * dt

    p_i.v .+= p_i.a * dt

    #reinitialize
    p_i.f.*= 0.

    return p_i
end
function overdamped_evolver!(p_i::Hexbug, t, dt)

    #evolve
    p_i.x .+= p_i.v * dt
    p_i.xuw.+= p_i.v * dt
    p_i.v .= p_i.f/p_i.zeta[1]

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end
function inertial_evolver!(p_i::PolarParticle3d, t, dt)


    #evolve
    p_i.a .= p_i.f/p_i.m[1]

    p_i.x .+= p_i.v * dt
    p_i.xuw.+= p_i.v * dt
    p_i.v .+= p_i.a * dt

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end

function overdamped_evolver!(p_i::PolarParticle3d, t, dt)

    #evolve
    p_i.x .+= p_i.v * dt
    p_i.xuw.+= p_i.v * dt
    p_i.v .= p_i.f/p_i.zeta[1]

    p_i.p .+= p_i.q * dt
    p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end

function overdamped_evolver!(p_i::ConfinedPolarParticle3d, t, dt)


    #evolve
    p_i.p .+= p_i.q * dt
    #p_i.p .=normalize(p_i.p)


    #reinitialize
    p_i.q.*= 0.
    p_i.f.*= 0.

    return p_i
end

function overdamped_evolver!(p_i::PolarParticle2d, t, dt)
    #Evolve
    p_i.x.+= p_i.v * dt
    
    p_i.xuw.+= p_i.v * dt

    p_i.v.= p_i.f/p_i.zeta[1]
    p_i.θ.+= p_i.ω * dt

    #reinitalize
    p_i.ω.*= 0.
    p_i.f.*= 0.

    return p_i

end



function overdamped_evolver!(p_i::PolarParticle2dN, t, dt)

    if p_i.n[1]>0
        p_i.f.+= p_i.fn/p_i.n[1]
    end

    #Evolve
    p_i.x.+= p_i.v * dt
    p_i.xuw.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta[1]
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
    p_i.xuw.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta[1]
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
    p_i.xuw.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta[1]
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
    p_i.xuw.+= p_i.v * dt
    p_i.v.= p_i.f/p_i.zeta[1]
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
function overdamped_evolver!(field::GeneralField2d, t, dt)

    field.C.+= field.Cv*dt

    field.Cv.= field.Cf

    #reinitalize
    field.Cf.*= 0.
    return field
end