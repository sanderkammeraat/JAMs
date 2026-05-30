
include("../analysis/AnalysisPipeline.jl")


base_folder = "/Users/kammeraat/surfdrive/ActivePolygonClusters/simulations/for_inference_v14_exp_rep/phi_0p01/simdata"
raw_data = load_file(joinpath(base_folder,"raw_data.h5"))

using DataFrames
using CSV


df = DataFrame(frame = Int[], particle=Int64[], angle_mapped=Float64[], x_um = Float64[], y_um=Float64[], vx_ums = Float64[], vy_ums=Float64[])


frame_numbers = 1:length(raw_data["frames"])


for i in frame_numbers


    frame = raw_data["frames"][string(i)]
    ids = frame["id"]


    for id in ids


        px = frame["px"][id]

        py = frame["py"][id]

        angle_id = angle(px+1im*py)
 
        x_um_id = frame["x"][id] * 5 #um

        y_um_id = frame["y"][id] * 5 #um

        vx_um_id = frame["vx"][id] * 5 #ums

        vy_um_id = frame["vy"][id] * 5 #ums

        #start at frame 0 and particle 0
        push!(df, [i-1, id-1 , angle_id, x_um_id, y_um_id, vx_um_id, vy_um_id])


    end
end

CSV.write(joinpath(base_folder, "traj.csv"), df)

close(raw_data)

