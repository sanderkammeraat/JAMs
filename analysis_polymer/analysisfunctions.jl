function save_data(data, path)

    save(joinpath(path,"analysis.jld2"), data)

end


function average_velocity(v_x, v_y, numb_frames, numb_particles, t_stop)
  
    v = zeros(numb_frames)
    dt = t_stop/numb_frames
    time = convert(Vector{Float64}, 0:dt:t_stop)


    for i in 1:numb_frames
        v[i] += sum(sqrt.(v_x[i].^2 + v_y[i].^2))/numb_particles
    end

    while length(v) != length(time)
        if length(v) > length(time)
            pop!(v)
        else
            pop!(time)
        end
    end

    return Dict("average_velocity" => v, "average_velocity_time" => time)
end



function MSD(x, y, numb_frames, numb_particles, t_stop)

    dt = t_stop/numb_frames
    MSD = zeros(numb_frames-1)
    time = convert(Vector{Float64}, 0:dt:t_stop-dt)

    for i in eachindex(MSD)
        MSD[i] += sum((x[i+1:numb_frames,:] - x[1:numb_frames-i,:]).^2 + (y[i+1:numb_frames,:] - y[1:numb_frames-i,:]).^2)/numb_particles/(numb_frames-i)
    end


    while length(MSD) != length(time)
        if length(MSD) > length(time)
            pop!(MSD)
        else
            pop!(time)
        end
    end

    return Dict("MSD" => MSD, "MSD_time" => time)
end



function basic_MSD(x, y, numb_frames, numb_particles, t_stop)
    
    dt = t_stop/numb_frames
    MSD = zeros(numb_frames)
    time = convert(Vector{Float64}, 0:dt:t_stop)
    x_0, y_0 = x[1,:], y[1,:]

    for i in eachindex(MSD)
        MSD[i] += sum((x[i,:] .- x_0).^2 + (y[i,:] .- y_0).^2)/numb_particles
    end

    while length(MSD) != length(time)
        if length(MSD) > length(time)
            pop!(MSD)
        else
            pop!(time)
        end
    end

    return Dict("basic_MSD" => MSD, "basic_MSD_time" => time)
end


function end_to_end_distance(x, y, pol_id, id_in_pol, numb_frames, Npol, N, t_stop)

    dt = t_stop/numb_frames
    time = convert(Vector{Float64}, 0:dt:t_stop)
    length_polymer = maximum(id_in_pol)
    end_to_end_distance = zeros(numb_frames)

    for i in 1:numb_frames
        for j in 1:Npol
            indices = findall(pol_id[pol_id=j])
            part_id = eachindex(id_in_pol[indices])
            x_begin = x[part_id]

            end_to_end_distance[i] += sqrt((x[(j-1)*length_polymer+1] - x[j*length_polymer])^2 + (y[(j-1)*length_polymer+1] - y[j*length_polymer])^2)/Npol
        end
    end

    while length(end_to_end_distance) != length(time)
        if length(end_to_end_distance) > length(time)
            pop!(end_to_end_distance)
        else
            pop!(time)
        end
    end

    return Dict("e_to_e_dist" => end_to_end_distance, "e_to_e_dist_time" => time)
end


function radius_of_gyration(x, y, pol_id, id_in_pol, numb_frames, numb_particles)

    length_polymer = maximum(id_in_pol)
    polymer_size = convert(Int64, length_polymer)
    numb_polymers = numb_particles/length_polymer
    R_2 = zeros(numb_frames)

    center_x = zeros(numb_frames, numb_particles)
    center_y = zeros(numb_frames, numb_particles)

    for i in 1:numb_frames

        for j in 1:numb_polymers
            center_x[i, convert(Int64, (j-1)*polymer_size+1):convert(Int64, j*polymer_size)] .+= sum(x[i, convert(Int64, (j-1)*polymer_size+1):convert(Int64, j*polymer_size)],dims=2)/length_polymer
            center_y[i, convert(Int64, (j-1)*polymer_size+1):convert(Int64, j*polymer_size)] .+= sum(x[i, convert(Int64, (j-1)*polymer_size+1):convert(Int64, j*polymer_size)],dims=2)/length_polymer
        end
        
    R_2 += sum((x .- center_x).^2 + (y .- center_y).^2, dims=2)/numb_particles
    end
    return R_2
end


