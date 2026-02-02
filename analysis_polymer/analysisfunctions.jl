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
    polymer_size = convert(Int64, N/Npol)
    end_to_end_distance = zeros(numb_frames)

    for i in 1:numb_frames
        for j in 1:Npol

            pol_indices = findall(index -> index == j, pol_id)

            begin_parts = findall(index -> index == 1, id_in_pol)
            end_parts = findall(index -> index == polymer_size, id_in_pol)

            begin_part_id = intersect(pol_indices, begin_parts)[1]
            end_part_id = intersect(pol_indices, end_parts)[1]

            end_to_end_distance[i] += sqrt((x[begin_part_id] - x[end_part_id])^2 + (y[begin_part_id] - y[end_part_id])^2)/polymer_size
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


function radius_of_gyration(x, y, pol_id, numb_frames, Npol, N, t_stop)

    dt = t_stop/numb_frames
    time = convert(Vector{Float64}, 0:dt:t_stop)
    polymer_size = convert(Int64, N/Npol)
    R_2 = zeros(numb_frames)


    for frame in 1:numb_frames

        for i in 1:Npol
            pol_indices = findall(index -> index == i, pol_id)
            for index_1 in pol_indices
                for index_2 in pol_indices
                    if index_1 != index_2
                        R_2[frame] += 1/(2*polymer_size^2) * ((x[index_1] - x[index_2]) .^2 .+ (y[index_1] - y[index_2]) .^2) / Npol
                    end
                end
            end
        end
    end
    
    while length(R_2) != length(time)
        if length(R_2) > length(time)
            pop!(R_2)
        else
            pop!(time)
        end
    end

    return Dict("R_2" => R_2, "R_2_time" => time)
end


