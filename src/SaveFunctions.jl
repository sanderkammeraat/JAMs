


function save_2d_polar!(file, current_particle_state, current_field_state, n, Tsave, t)

    preamble = string(n)*"/"


    file[preamble*"n"] = n

    file[preamble*"t"] = t

    file[preamble*"id"] = [p_i.id for p_i in current_particle_state]

    file[preamble*"type"] = [p_i.type for p_i in current_particle_state]

    file[preamble*"structtype"] = [string(nameof(typeof(p_i))) for p_i in current_particle_state]

    file[preamble*"x"] = [p_i.x[1] for p_i in current_particle_state]
    file[preamble*"y"] = [p_i.x[2] for p_i in current_particle_state]

    file[preamble*"vx"] = [p_i.v[1] for p_i in current_particle_state]
    file[preamble*"vy"] = [p_i.v[2] for p_i in current_particle_state]

    file[preamble*"px"] = [p_i.p[1] for p_i in current_particle_state]
    file[preamble*"py"] = [p_i.p[2] for p_i in current_particle_state]

    file[preamble*"qx"] = [p_i.q[1] for p_i in current_particle_state]
    file[preamble*"qy"] = [p_i.q[2] for p_i in current_particle_state]

    return file


end