
#new, local_evolver and global_evolver

struct overdamped_pairdis_evolver
    η::Float64
    rfact::Float64 #should be smaller than the cut off radius
end


struct polymer_p_set_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_xvf_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_pq_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_pq_xyc_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_2d_shape_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_θω_evolver
    ontypes::Union{Int64,Vector{Int64}}
end

struct overdamped_deformable_ellipse_evolver
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
#Note that a potential Itô relaxation term is taken care of by renormalizing p_i.p
function evolve_locally!(p_i, t, dt, dofevolver::overdamped_pq_evolver)

    if p_i.type[1] in dofevolver.ontypes

        #evolve
        p_i.p .+= cross(p_i.q,p_i.p) * dt
        p_i.p .=normalize(p_i.p)

        #reinitialize
        p_i.q.*= 0.
    end
    return p_i
end
#constrained to be in xy
function evolve_locally!(p_i, t, dt, dofevolver::overdamped_pq_xyc_evolver)

    if p_i.type[1] in dofevolver.ontypes


        dθ = p_i.q[3] * dt
        #evolve
        pxc = copy(p_i.p[1])
        pyc = copy(p_i.p[2])
        p_i.p[1] = cos(dθ) * pxc  - sin(dθ) *  pyc
        p_i.p[2]= sin(dθ) * pxc  + cos(dθ) *  pyc
        p_i.p .=normalize(p_i.p)

        #reinitialize
        p_i.q.*= 0.
    end
    return p_i
end

function evolve_locally!(p_i, t, dt, dofevolver::overdamped_2d_shape_evolver)

    if p_i.type[1] in dofevolver.ontypes

        for j = 1:size(p_i.xe)[1]

             xen =  p_i.xo[j,1]* p_i.p[1] -   p_i.xo[j,2] * p_i.p[2] 
             yen = p_i.xo[j,1]* p_i.p[2]  +  p_i.xo[j,2] * p_i.p[1]

             p_i.xe[j,1] = p_i.x[1] + xen
             p_i.xe[j,2] = p_i.x[2] + yen
        end
    end
    return p_i
end



function evolve_locally!(p_i, t, dt, dofevolver::overdamped_θω_evolver)

    if p_i.type[1] in dofevolver.ontypes
    #evolve
    p_i.θ.+= p_i.ω .* dt

    #reinitalize
    p_i.ω.*= 0.
    end

    return p_i
end

function evolve_locally!(p_i, t, dt, dofevolver::overdamped_deformable_ellipse_evolver)
    if p_i.type[1] in dofevolver.ontypes

        stress = p_i.stress
        R_0 = p_i.sigma_0[1]
        Lmda_0 = p_i.Lambda0
        Lmda = p_i.Lambda
        tau = p_i.tau[1]
        mu = p_i.mu[1]
        K = p_i.K[1]

        dLmda_dt = @MMatrix [0.0 0.0; 0.0 0.0]
        
        dLmda_dt[1,1] = -1/tau*(Lmda[1,1]-Lmda_0[1,1]) + R_0/(4*tau*mu*K)*((mu+K)*stress[1,1] + (mu-K)*stress[2,2])
        dLmda_dt[2,2] = -1/tau*(Lmda[2,2]-Lmda_0[2,2]) + R_0/(4*tau*mu*K)*((mu+K)*stress[2,2] + (mu-K)*stress[1,1])
        dLmda_dt[1,2] = -1/tau*(Lmda[1,2]-Lmda_0[1,2]) + R_0/(2*tau*mu)*stress[1,2]
        dLmda_dt[2,1] = -1/tau*(Lmda[2,1]-Lmda_0[2,1]) + R_0/(2*tau*mu)*stress[2,1]

        p_i.Lambda.+= dLmda_dt .*dt

        lmda_major = 0.5*(Lmda[1,1]+Lmda[2,2]) + sqrt(0.25*(Lmda[1,1]-Lmda[2,2])^2 + Lmda[1,2]^2)
        lmda_minor = 0.5*(Lmda[1,1]+Lmda[2,2]) - sqrt(0.25*(Lmda[1,1]-Lmda[2,2])^2 + Lmda[1,2]^2)

        theta_def = 0.5*atan(2*Lmda[1,2],(Lmda[1,1]-Lmda[2,2])) 
        #theta_def = theta_def - np.pi*round((theta_def-p_i.θ)/pi) 

        #Put in minimal angle 
        theta_current = angle.(p_i.p[1] + 1im * p_i.p[2])

        dθ = p_i.q[3] * dt + (theta_def - theta_current)
        #evolve
        pxc = copy(p_i.p[1])
        pyc = copy(p_i.p[2])
        p_i.p[1] = cos(dθ) * pxc  - sin(dθ) *  pyc
        p_i.p[2]= sin(dθ) * pxc  + cos(dθ) *  pyc
        p_i.p .=normalize(p_i.p)

        #reinitialize
        p_i.q.*= 0.
        
        sn = p_i.p[2]
        cs = p_i.p[1]

        p_i.Lambda[1,1] = cs^2*lmda_major+sn^2*lmda_minor
        p_i.Lambda[2,2] = sn^2*lmda_major+cs^2*lmda_minor
        p_i.Lambda[1,2] = sn*cs*(lmda_major-lmda_minor)
        p_i.Lambda[2,1] = p_i.Lambda[1,2]
    
        p_i.stress.*=0
    end
    return p_i
end


function evolve_globally!(current_particle_state, current_field_state, system, cells, stencils, dt, dofevolver::overdamped_pairdis_evolver)

    dims = length(system.sizes)
    Np = length(current_particle_state)
    id = SMatrix{dims, dims}(1I)
    id_global = SMatrix{dims, dims}(dofevolver.η*I) 


    ax1_ind = Int64[]
    ax2_ind = Int64[]
    vals = Float64[]

    self = @MMatrix zeros(Float64,dims, dims)

    #We are going to solve Mx=b with b the forces (excluding pair-dissipation) and M = I \zeta + M_pairdiss.
    # For that we need to construct M, first preallocate (make sparse?)
    for i in eachindex(current_particle_state)

        #self part of matrix M
        
        p_i = current_particle_state[i]

        dx = @MVector zeros(Float64,length(p_i.x))
        neighbours= get_neighbours(p_i,cells,stencils)
        if !isnothing(neighbours)
    
            for n in neighbours
    
                if i!=n
                    p_j = current_particle_state[n]
    
                    dx = minimal_image_difference!(dx, p_i.x, p_j.x, system.sizes, system.Periodic)
    
                    dxn = norm(dx)

                    d2a = p_i.R[1]+p_j.R[1]

                    range = dofevolver.rfact*d2a::Float64
                    
                    if dxn<=range #In contact!


                        for (im, Mi) in pairs(dims*(i-1)+1:dims*i)

                            for (jm, Mj) in pairs(dims*(n-1)+1:dims*n)

                                if im==jm
                                    push!(ax1_ind,Mi)
                                    push!(ax2_ind,Mj)
                                    push!(vals, -id_global[im,jm])
                                    self[im,jm]+=id_global[im,jm]
                            
                                end
                            end
                        end
                        #M[dims*(i-1)+1:dims*i, dims*(n-1)+1:dims*n].-=id_global

                        #M[dims*(i-1)+1:dims*i, dims*(i-1)+1:dims*i].+= id_global
                    end
                end
                if i==n
                    for (im, Mi) in pairs(dims*(i-1)+1:dims*i)
                        self[im,im]+= id[im,im]*p_i.zeta[1]
                    end
                end

            end    
            for (im, Mi) in pairs(dims*(i-1)+1:dims*i)
                push!(ax1_ind,Mi)
                push!(ax2_ind,Mi)
                push!(vals, self[im,im])
            end
            #reinitialize
            self.*=0.

            
        end
    end

    M = sparse(ax1_ind, ax2_ind, vals, dims*Np, dims*Np)

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

function evolve_globally!(current_particle_state, current_field_state, system, cells, stencils, dt, dofevolver::polymer_p_set_evolver)

    Threads.@threads for i in eachindex(current_particle_state)

        p_i = current_particle_state[i]

        #discard current polarity vector
        p_i.p.*=0

        dx = @MVector zeros(Float64,length(p_i.x))
        neighbours= get_neighbours(p_i,cells,stencils)
        if !isnothing(neighbours)
    
            for n in neighbours
    
                if i!=n
                    p_j = current_particle_state[n]
    
                    dx = minimal_image_difference!(dx, p_i.x, p_j.x, system.sizes, system.Periodic)
    
                    dxn = norm(dx)
                    
                    if p_i.type[1] in dofevolver.ontypes && p_j.type[1] in dofevolver.ontypes

                            if p_i.pol_id[1]==p_j.pol_id[1]

                                if p_j.id_in_pol[1]==p_i.id_in_pol[1]+1

                                    p_i.p.+=dx/dxn

                                elseif  p_j.id_in_pol[1]==p_i.id_in_pol[1]-1

                                    p_i.p.+= -dx/dxn
                                end
                            end
                        end
                end
            end    
        end
        p_i.p.= normalize(p_i.p)
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



#old, no longer functional, focus on particle

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

struct overdamped_CCvCf_evolver
    ontypes::Union{Int64,Vector{Int64}}
end
    
function evolve_field!(field, t, dt, dofevolver::overdamped_CCvCf_evolver)

    if field.type in dofevolver.ontypes
        field.C.+= field.Cv*dt

        field.Cv.= field.Cf

        #reinitalize
        field.Cf.*= 0.
    end
    return field
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