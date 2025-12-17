


function average_velocity(v_x, v_y, numb_frames, numb_particles, sliding_window)
  
    v = zeros(numb_frames-sliding_window+1)

    for j in 1:sliding_window
        v[1] += sum(sqrt.(v_x[j].^2 + v_y[j].^2))/numb_particles
    end

    for i in 2:numb_frames-sliding_window+1
        v[i] += v[i-1] + sum(sqrt.(v_x[i+sliding_window-1].^2 + v_y[i+sliding_window-1].^2))/numb_particles - sum(sqrt.(v_x[i-1].^2 + v_y[i-1].^2))/numb_particles
    end

    return v/sliding_window
end



function MSD(x, y, numb_frames, numb_particles)

    MSD = zeros(convert(Int64, floor((numb_frames-1)/4)))

    for Δt in 1:4:numb_frames-17
        MSD[convert(Int64, floor((Δt-1)/4))+1] += sum((x[Δt+1:numb_frames,:] - x[1:numb_frames-Δt,:]).^2 + (y[Δt+1:numb_frames,:] - y[1:numb_frames-Δt,:]).^2)/numb_particles/(numb_frames-Δt)
    end

    return MSD

end

function plot_MSD!(data, end_time)

    f = Figure()
    ax = Axis(f[1, 1], xscale=log10, yscale=log10)
    time = 1:length(data)
    ylims!(ax, 1e-3,maximum(data)*1.5)


    lines!(ax, time*end_time/length(data), data)
    display(f)
end


function end_to_end_distance(x, y, pol_id, id_in_pol, numb_frames, numb_particles)

    length_polymer = max(id_in_pol)
    numb_polymers = numb_particles/
    end_to_end_distance = zeros(numb_frames)
    for i in 1:numb_polymers
        end_to_end_distance[i] = 2
    end
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

function plot_radius_of_gyration!(multiple, datas, p, end_time)

    if !multiple

        f = Figure()

        ax = Axis(f[1, 1], xscale=log10, yscale=log10)
        time = 1:length(datas)
        lines!(time*end_time/length(data), datas)
        display(f)

    
    else
        f = Figure()
        ax = Axis(f[1, 1])

        i=1
        for data in datas
            activity = p[i]
            time = 1:length(data)
            lines!(time*end_time/length(data), data, label="p=$activity")
            i+=1
        end
    end

    Legend(f[1,2], ax, "activity")
    display(f)

end