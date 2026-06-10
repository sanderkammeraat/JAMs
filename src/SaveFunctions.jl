#Refer to https://juliaio.github.io/HDF5.jl/stable/#Opening-and-closing-files
# for allowed datastructures


#Currently only supporting a single field only
function save_single_2d_field!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t,framecounter)



    current_frame_group["field_n"] = n

    current_frame_group["field_t"] = t


    current_frame_group["field_id"] = [field.id for field in current_field_state]

    current_frame_group["field_type"] = [field.type for field in current_field_state]
    current_frame_group["field_bin_centers_x"] = current_field_state[1].bin_centers[1]
    current_frame_group["field_bin_centers_y"] = current_field_state[1].bin_centers[2]

    current_frame_group["field_C"] = current_field_state[1].C


    return current_frame_group

end

function save_single_2d_field_10!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t,framecounter)


    if (n-1)%(10*Tsave)==0
        current_frame_group["field_n"] = n

        current_frame_group["field_t"] = t


        current_frame_group["field_id"] = [field.id for field in current_field_state]

        current_frame_group["field_type"] = [field.type for field in current_field_state]
        current_frame_group["field_bin_centers_x"] = current_field_state[1].bin_centers[1]
        current_frame_group["field_bin_centers_y"] = current_field_state[1].bin_centers[2]

        current_frame_group["field_C"] = current_field_state[1].C
    end

    return current_frame_group

end

function save_2d_polar_p!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t,framecounter)




    current_frame_group["n"] = n

    current_frame_group["t"] = t

    current_frame_group["id"] = [p_i.id[1] for p_i in current_particle_state]

    current_frame_group["type"] = [p_i.type[1] for p_i in current_particle_state]

    current_frame_group["v0"] = [p_i.v0[1] for p_i in current_particle_state]

    current_frame_group["Dr"] = [p_i.Dr[1] for p_i in current_particle_state]

    current_frame_group["R"] = [p_i.R[1] for p_i in current_particle_state]

    current_frame_group["structtype"] = [string(nameof(typeof(p_i))) for p_i in current_particle_state]

    current_frame_group["x"] = [p_i.x[1] for p_i in current_particle_state]
    current_frame_group["y"] = [p_i.x[2] for p_i in current_particle_state]

    current_frame_group["xuw"] = [p_i.xuw[1] for p_i in current_particle_state]
    current_frame_group["yuw"] = [p_i.xuw[2] for p_i in current_particle_state]

    current_frame_group["vx"] = [p_i.v[1] for p_i in current_particle_state]
    current_frame_group["vy"] = [p_i.v[2] for p_i in current_particle_state]

    current_frame_group["px"] = [p_i.p[1] for p_i in current_particle_state]
    current_frame_group["py"] = [p_i.p[2] for p_i in current_particle_state]

    current_frame_group["qx"] = [p_i.q[1] for p_i in current_particle_state]
    current_frame_group["qy"] = [p_i.q[2] for p_i in current_particle_state]

    return current_frame_group


end




function save_2d_shape_polar_p!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t,framecounter)

    current_frame_group["n"] = n

    current_frame_group["t"] = t

    current_frame_group["id"] = [p_i.id[1] for p_i in current_particle_state]

    current_frame_group["type"] = [p_i.type[1] for p_i in current_particle_state]

    current_frame_group["v0"] = [p_i.v0[1] for p_i in current_particle_state]

    current_frame_group["Dr"] = [p_i.Dr[1] for p_i in current_particle_state]

    current_frame_group["R"] = [p_i.R[1] for p_i in current_particle_state]

    current_frame_group["xe_x"] = [p_i.xe[j,1] for p_i in current_particle_state for j=1:size(p_i.xe)[1] ]

    current_frame_group["xe_y"] = [p_i.xe[j,2] for p_i in current_particle_state for j=1:size(p_i.xe)[1] ]

    current_frame_group["structtype"] = [string(nameof(typeof(p_i))) for p_i in current_particle_state]

    current_frame_group["x"] = [p_i.x[1] for p_i in current_particle_state]
    current_frame_group["y"] = [p_i.x[2] for p_i in current_particle_state]

    current_frame_group["xuw"] = [p_i.xuw[1] for p_i in current_particle_state]
    current_frame_group["yuw"] = [p_i.xuw[2] for p_i in current_particle_state]

    current_frame_group["vx"] = [p_i.v[1] for p_i in current_particle_state]
    current_frame_group["vy"] = [p_i.v[2] for p_i in current_particle_state]

    current_frame_group["px"] = [p_i.p[1] for p_i in current_particle_state]
    current_frame_group["py"] = [p_i.p[2] for p_i in current_particle_state]

    current_frame_group["qx"] = [p_i.q[1] for p_i in current_particle_state]
    current_frame_group["qy"] = [p_i.q[2] for p_i in current_particle_state]

    return current_frame_group


end
function save_2d_polymer_polar_p!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t,framecounter)




    current_frame_group["n"] = n

    current_frame_group["t"] = t

    current_frame_group["id"] = [p_i.id[1] for p_i in current_particle_state]

    current_frame_group["pol_id"] = [p_i.pol_id[1] for p_i in current_particle_state]

    current_frame_group["id_in_pol"] = [p_i.id_in_pol[1] for p_i in current_particle_state]

    current_frame_group["pol_N"] = [p_i.pol_N[1] for p_i in current_particle_state]

    current_frame_group["type"] = [p_i.type[1] for p_i in current_particle_state]

    current_frame_group["v0"] = [p_i.v0[1] for p_i in current_particle_state]

    current_frame_group["Dr"] = [p_i.Dr[1] for p_i in current_particle_state]

    current_frame_group["R"] = [p_i.R[1] for p_i in current_particle_state]

    current_frame_group["structtype"] = [string(nameof(typeof(p_i))) for p_i in current_particle_state]

    current_frame_group["x"] = [p_i.x[1] for p_i in current_particle_state]
    current_frame_group["y"] = [p_i.x[2] for p_i in current_particle_state]

    current_frame_group["xuw"] = [p_i.xuw[1] for p_i in current_particle_state]
    current_frame_group["yuw"] = [p_i.xuw[2] for p_i in current_particle_state]

    current_frame_group["vx"] = [p_i.v[1] for p_i in current_particle_state]
    current_frame_group["vy"] = [p_i.v[2] for p_i in current_particle_state]

    current_frame_group["px"] = [p_i.p[1] for p_i in current_particle_state]
    current_frame_group["py"] = [p_i.p[2] for p_i in current_particle_state]

    current_frame_group["qx"] = [p_i.q[1] for p_i in current_particle_state]
    current_frame_group["qy"] = [p_i.q[2] for p_i in current_particle_state]

    return current_frame_group


end

function save_2d_polar_θ!(current_frame_group, current_particle_state, current_field_state, n, Tsave, t, framecounter)


    current_frame_group["n"] = n

    current_frame_group["t"] = t

    current_frame_group["id"] = [p_i.id[1] for p_i in current_particle_state]

    current_frame_group["type"] = [p_i.type[1] for p_i in current_particle_state]

    current_frame_group["v0"] = [p_i.v0[1] for p_i in current_particle_state]

    current_frame_group["Dr"] = [p_i.Dr[1] for p_i in current_particle_state]

    current_frame_group["R"] = [p_i.R[1] for p_i in current_particle_state]

    current_frame_group["structtype"] = [string(nameof(typeof(p_i))) for p_i in current_particle_state]

    current_frame_group["x"] = [p_i.x[1] for p_i in current_particle_state]
    current_frame_group["y"] = [p_i.x[2] for p_i in current_particle_state]

    current_frame_group["xuw"] = [p_i.xuw[1] for p_i in current_particle_state]
    current_frame_group["yuw"] = [p_i.xuw[2] for p_i in current_particle_state]

    current_frame_group["vx"] = [p_i.v[1] for p_i in current_particle_state]
    current_frame_group["vy"] = [p_i.v[2] for p_i in current_particle_state]

    current_frame_group["θ"] = [p_i.θ[1] for p_i in current_particle_state]
    current_frame_group["ω"] = [p_i.ω[1] for p_i in current_particle_state]

    return current_frame_group


end