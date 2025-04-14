

#Multi-particle forces: Neighbourhood
struct swarm_pos_force <: NeighbourhoodForce
    ontypes::Union{Int64,Vector{Int64}}
    J::Float64
end

struct swarm_angular_force<:NeighbourhoodForce
    ontypes::Union{Int64,Vector{Int64}}
    K::Float64
end

struct Vicsek_align_force<:NeighbourhoodForce
    ontypes::Union{Int64,Vector{Int64}}
    r::Float64

end


#Template by defining a contribute function for general type

function contribute_neighbourhood_force!(p_i, p_js, dxs, dxns, t, dt, force,rngs_particles)

    #p_js is an array containing the particles that are within rcut_pair_global

    #dxs is an array containing the dx for each p_j in p_js

    #dxns is an array containing the dxn for each p_j in p_js

    if p_i.type[1] in force.ontypes
        neighbourhood = [j for j in eachindex(p_js) if (p_js[j].type[1] in force.ontypes) ]    
        N = lenght(neighbourhood)    
        #check against something else to create smaller subset and calculate
    end
    return p_i
end


function contribute_neighbourhood_force!(p_i, p_js, dxs, dxns, t, dt, force::swarm_pos_force,rngs_particles)
    if p_i.type[1] in force.ontypes
        neighbourhood = [j for j in eachindex(p_js) if (p_js[j].type[1] in force.ontypes) ]  
        N = length(neighbourhood)
        for j in neighbourhood
            @views p_i.f.+= 1/N * (dxs[j]/dxns[j] * (1 + force.J*cos(p_js[j].ϕ[1]-p_i.ϕ[1]) ) - dxs[j]/dxns[j]^2)  
        end
    end 
    return p_i
end

function contribute_neighbourhood_force!(p_i, p_js, dxs, dxns, t, dt,  force::swarm_angular_force,rngs_particles)
    if p_i.type[1] in force.ontypes

        neighbourhood = [j for j in eachindex(p_js) if (p_js[j].type[1] in force.ontypes) ]  
        N = length(neighbourhood)

        for j in neighbourhood
            p_i.ψ.+= 1/N * force.K * sin(p_js[j].ϕ[1]-p_i.ϕ[1])/dxns[j]
        end
    end
    return p_i

end



function contribute_neighbourhood_force!(p_i, p_js, dxs, dxns, t, dt,  force::Vicsek_align_force,rngs_particles)
    if p_i.type[1] in force.ontypes
        neighbourhood = [j for j in eachindex(p_js) if ((p_js[j].type[1] in force.ontypes) && (dxns[j]< force.r)) ] 
        N = length(neighbourhood) 
        for j in neighbourhood
            println(N)
            p_i.ω[1] += p_js[j].θ[1]/dt/N
        end
    end
    return p_i
end



