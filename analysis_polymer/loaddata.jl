using JLD2


function load_file(file_location)

    file = jldopen(file_location, "r", iotype=IOStream)

end


struct simulation_data
    x
    y
    vx
    vy
    px
    py
    id::Vector{Int64}
    pol_id::Vector{Int64}
    id_in_pol::Vector{Int64}
    R::Vector{Float64}
    numb_frames::Int64
    dt::Float64
    Tsave::Int64
    ksd::Float64
    kbend::Float64
    kstretch::Float64
    fstretch::Float64
    #p::Float64
    #kperp::Float64
    #kpar::Float64
    Npol::Int64
    N::Int64
end



function get_data(file_location)

    file = load_file(file_location)
    

    numb_frames = length(file["frames"])
    N = length(file["frames/1/xuw"])
    pol_id = file["frames/1/pol_id"]
    id_in_pol = file["frames/1/id_in_pol"]


    v_x = zeros(numb_frames, N)
    v_y = zeros(numb_frames, N)

    x = zeros(numb_frames, N)
    y = zeros(numb_frames, N)

    p_x = zeros(numb_frames, N)
    p_y = zeros(numb_frames, N)

    for i in 1:numb_frames

        x[i,:] = file["frames/$i/xuw"]
        y[i,:] = file["frames/$i/yuw"]
        v_x[i,:] = file["frames/$i/vx"]
        v_y[i,:] = file["frames/$i/vy"]
        p_x[i,:] = file["frames/$i/px"]
        p_y[i,:] = file["frames/$i/py"]

    end

    id = file["frames/5/id"]
    pol_id = file["frames/5/pol_id"]
    id_in_pol = file["frames/5/id_in_pol"]
    R = file["frames/5/R"]
    dt = file["integration_info/dt"]
    Tsave = file["integration_info/Tsave"]
    ksd = file["system/forces/pair/polymer_exterior_soft_disk_force/karray"]
    kbend = file["system/forces/pair/polymer_harmonic_bend_force/karray"]
    kstretch = file["system/forces/pair/polymer_harmonic_stretch_force/karray"]
    fstretch = file["system/forces/pair/polymer_harmonic_stretch_force/farray"]
    #p = file["system/forces/pair/polymer_pairAN_force/parray"]
    #kperp = file["system/forces/pair/polymer_pairAN_force/k_perp"]
    #kpar = file["system/forces/pair/polymer_pairAN_force/k_par"]
    Npol = maximum(pol_id)

    close(file)

    return simulation_data(x, y, v_x, v_y, p_x, p_y, id, pol_id, id_in_pol, R, numb_frames, dt, Tsave, ksd, kbend, kstretch, fstretch, Npol, N)
end