struct field_propulsion_force<:FieldForce
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    v0offset::Float64
    
end
struct field_propulsion_3d_force<:FieldForce
    ontypes::Union{Int64,Vector{Int64}}
    consumption::Float64
    v0offset::Float64
end

#Template. Note that field forces are not multithreaded and can therefore change the state of the field as well, without worrying about threadsafety.
function contribute_field_force!(p_i,field_j,field_indices, t, dt, force,rngs_particles)
    if p_i.type[1] in force.ontypes && field_j.type in force.ontypes
        #Do calculation
    end

    return p_i, field_j

end


function contribute_field_force!(p_i,field_j,field_indices, t, dt, force::field_propulsion_force,rngs_particles)
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

function contribute_field_force!(p_i,field_j,field_indices, t, dt, force::field_propulsion_3d_force, rngs_particles)
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