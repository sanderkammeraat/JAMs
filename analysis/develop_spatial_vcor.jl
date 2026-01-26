



file_path = "/Users/kammeraat/mounting/data2_kammeraat/sa/statistics/hex_disordered/phi_1.3/Nlin_20/simdata/v0_0.01/Dr_0.01/J_0.1/seed_1/sa_raw_data.h5"


include("AnalysisPipeline.jl")


file = load_file(file_path)


frames = file["frames"]



t = file["integration_info"]["save_tax"]


Nint = sum(frames["1"]["type"].==1)
Nt = length(t)

Nt=2000

x = zeros(Nint, Nt)
y = zeros(Nint, Nt)

vx = zeros(Nint, Nt)
vy = zeros(Nint, Nt)



@views for i in 1:Nt
    x[:,i] .= extract_frame_data_for_type("x", 1, frames[string(i)])
    y[:,i] .= extract_frame_data_for_type("y", 1, frames[string(i)])
    vx[:,i] .= extract_frame_data_for_type("vx", 1, frames[string(i)])
    vy[:,i] .= extract_frame_data_for_type("vy", 1, frames[string(i)])
end


function spatial_v_correlation(binsize, maxbin_center, x,y, vx, vy)


    rbin_edges = prepend!(collect(range(start=0, step=binsize, stop=maxbin_center)),[0])

    rbin_edges2 = rbin_edges.^2

    rbin_centers = (rbin_edges[2:end] + rbin_edges[1:end-1])/2

    Nbin = length(rbin_centers)

    Nt = size(vx)[2]

    C = ones(Nt, Nbin)*NaN

    counts = zeros(Nt, Nbin)

     @showprogress dt = 1 desc="spatial correlation" showspeed=true Threads.@threads for i in 1:Nt

        vrms_2_i= mean( vx[:, i].^2 .+ vy[:, i].^2 )
        for v1 in 1:size(vx)[1]

            for v2 in v1:size(vx)[1]

                Δr2 = (x[v1, i]- x[v2,i])^2 +   (y[v1, i]- y[v2,i])^2

                for bin in eachindex(rbin_centers)
                    #inbounds, compiler may not know that rbin_centers is sh

                    @inbounds if (Δr2<= rbin_edges2[bin+1] && Δr2> rbin_edges2[bin]) || (v1==v2 && bin==1)

                        if isnan(C[i, bin])
                            C[i, bin]=0
                            counts[i, bin]=0
                        end

                        C[i, bin] +=( vx[v1, i]* vx[v2, i] + vy[v1, i] * vy[v2, i])/vrms_2_i
                        counts[i,bin] += 1
                    end

                end

            end
        end
    end

    C.=C./counts
    Cavg = mean(C, dims=1)[1,:]
    return Dict("rbc"=>rbin_centers,"rbe"=>rbin_edges,"Cavg"=> Cavg)

end


min_t_ind = 1000
spatial_vcor = spatial_v_correlation(3, 100, x[:,min_t_ind:end], y[:,min_t_ind:end], vx[:,min_t_ind:end], vy[:,min_t_ind:end])








begin 
using GLMakie
GLMakie.activate!()

f = Figure()

ax = Axis(f[1,1])

scatter!(ax, spatial_vcor["rbc"], spatial_vcor["Cavg"])

display(f)
end