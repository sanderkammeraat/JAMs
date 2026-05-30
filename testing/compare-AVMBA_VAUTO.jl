

using Glob
using ProgressBars
using CSV
using DataFrames

function parse_dat(load_filepath)

    full_keys = CSV.read(load_filepath, limit=1, delim=" ", ignorerepeated=true, DataFrame, header=false)

    #Strip off 'keys:' from first line of output. 
    keys =Vector{String}(full_keys[1,2:end])

    df = CSV.read(load_filepath, skipto=2, delim=" ", ignorerepeated=true, DataFrame, header=keys)

    return df
 
end

function scan_dat(load_folderpath, pattern)

    filepaths = glob(pattern, load_folderpath)

    return filepaths

end

function read_frames(load_filepaths)

    frames_data = []
    frames_name = []

    for (i, filepath) in ProgressBar(pairs(load_filepaths))

        df_i = parse_dat(filepath)
        push!(frames_data, df_i)
        
        frame_name_i = split(filepath,'/')[end]
        push!(frames_name, frame_name_i)

    end

    return frames_name, frames_data
end


fp_sim = scan_dat("/Volumes/T7_GREY/SAVM/v0_50/runs/v0_50/Jv0_1e-05/seed_0/autogen_AVM_simulation_step/","cell*.dat")



frames_name, frames_data = read_frames(fp_sim)


Nint = sum(frames_data[1].type.==1)
Nt =  length(frames_data)

x = zeros(Nint, Nt)
y = zeros(Nint, Nt)

vx = zeros(Nint, Nt)
vy = zeros(Nint, Nt)

#R = collect(frames_data[1].radius)
type = collect(frames_data[1].type)
@views for i in 1:Nt
    x[:,i] .=  frames_data[i].x[frames_data[i].type.==1]
    y[:,i] .=  frames_data[i].y[frames_data[i].type.==1]
    
    vx[:,i] .=  frames_data[i].vx[frames_data[i].type.==1]
    vy[:,i] .=  frames_data[i].vy[frames_data[i].type.==1]
end




include("../analysis/AnalysisPipeline.jl") 
include("../analysis/AnalysisFunctions.jl") 

t = collect(range(0,Nt-1))
VAUTO = auto_correlation(t, vx, vy, normalized=true)
begin
f = Figure()
ax = Axis(f[1,1])
scatter!(ax,VAUTO["t"]*0.1, VAUTO["Cavg"])
xlims!(ax,0,5)
display(f)
end

#PASS