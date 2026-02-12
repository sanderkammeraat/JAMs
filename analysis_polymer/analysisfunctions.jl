function average_velocity(v_x, v_y, numb_frames, numb_particles, t_stop)

    dt = t_stop/(numb_frames-1)
    time = convert(Vector{Float64}, 0:dt:t_stop)

    v = vec(sum(sqrt.(v_x.^2 .+ v_y.^2), dims=2) ./ numb_particles)

    return Dict("average_velocity" => v, "average_velocity_time" => time)
end



function MSD(x, y, numb_frames, numb_particles, t_stop)

    dt = t_stop/(numb_frames-1)
    MSD = zeros(numb_frames-1)
    time = convert(Vector{Float64}, 0:dt:t_stop-dt)

    for Δt in eachindex(MSD)
        for t in 1:20:numb_frames-Δt

            MSD[Δt] += sum((x[t+Δt:numb_frames] - x[t:numb_frames-Δt]).^2 .+ (y[t+Δt:numb_frames] - y[t:numb_frames-Δt]).^2)/numb_particles/length(1:20:numb_frames-Δt)
       
        end
    end

    return Dict("MSD" => MSD, "MSD_time" => time)
end


function polymer_MSD(x, y, pol_id, numb_frames, Npol, numb_particles, t_stop)

    dt = t_stop/(numb_frames-1)
    polymer_MSD = zeros(numb_frames-1)
    time = convert(Vector{Float64}, 0:dt:t_stop-dt)
    polymer_size = convert(Int64, numb_particles/Npol)


    polymer_x = zeros(numb_frames, Npol)
    polymer_y = zeros(numb_frames, Npol)

    for i in 1:Npol
        indices = findall(index -> index == i, pol_id)

        polymer_x[:, i] = sum(x[:, indices], dims=2) ./ polymer_size
        polymer_y[:, i] = sum(y[:, indices], dims=2) ./ polymer_size

    end

    for Δt in eachindex(1:numb_frames-1)
        for t in 1:20:numb_frames-Δt

            polymer_MSD[Δt] += sum((polymer_x[t+Δt:numb_frames, :] .- polymer_x[t:numb_frames-Δt]).^2 .+ (polymer_y[t+Δt:numb_frames] .- polymer_y[t:numb_frames-Δt]).^2)/Npol/length(1:10:numb_frames-Δt)
       
        end
    end

    return Dict("polymer_MSD" => polymer_MSD, "polymer_MSD_time" => time)
end



function basic_MSD(x, y, numb_frames, numb_particles, t_stop)
    
    dt = t_stop/(numb_frames-1)
    time = convert(Vector{Float64}, 0:dt:t_stop)
    x_0, y_0 = x[1,:], y[1,:]

    MSD = vec(sum((transpose(transpose(x) .- x_0).^2 .+ (transpose(y) .- y_0).^2), dims=2) ./ numb_particles)

    return Dict("basic_MSD" => MSD, "basic_MSD_time" => time)
end


function end_to_end_distance(x, y, pol_id, id_in_pol, numb_frames, Npol, N, t_stop)

    dt = t_stop/(numb_frames-1)
    time = convert(Vector{Float64}, 0:dt:t_stop)
    polymer_size = convert(Int64, N/Npol)
    end_to_end_distance = zeros(numb_frames)

    for j in 1:Npol

        pol_indices = findall(index -> index == j, pol_id)

        begin_parts = findall(index -> index == 1, id_in_pol)
        end_parts = findall(index -> index == polymer_size, id_in_pol)

        begin_part_id = intersect(pol_indices, begin_parts)[1]
        end_part_id = intersect(pol_indices, end_parts)[1]

        end_to_end_distance += vec(sqrt.((x[:, begin_part_id] .- x[:, end_part_id]).^2 .+ (y[:, begin_part_id] .- y[:, end_part_id]).^2) ./ polymer_size)
    end

    return Dict("e_to_e_dist" => end_to_end_distance, "e_to_e_dist_time" => time)
end


function radius_of_gyration(x, y, pol_id, numb_frames, Npol, N, t_stop)

    dt = t_stop/(numb_frames-1)
    time = convert(Vector{Float64}, 0:dt:t_stop)
    polymer_size = convert(Int64, N/Npol)
    R_2 = zeros(numb_frames)

    polymer_x = zeros(numb_frames, Npol)
    polymer_y = zeros(numb_frames, Npol)

    for i in 1:Npol
        indices = findall(index -> index == i, pol_id)

        polymer_x[:, i] = sum(x[:, indices], dims=2) ./ polymer_size
        polymer_y[:, i] = sum(y[:, indices], dims=2) ./ polymer_size

    end


    for i in 1:Npol
        pol_indices = findall(index -> index == i, pol_id)
        R_2 += vec(sum(((x[:, pol_indices] .- polymer_x[:, i]) .^2 .+ (y[:, pol_indices] .- polymer_y[:, i]) .^2), dims=2) ./ N)
    end

    return Dict("R_2" => R_2, "R_2_time" => time)
end


