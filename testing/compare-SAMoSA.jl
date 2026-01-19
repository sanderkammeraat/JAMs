

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


fp_sim = scan_dat("/Users/kammeraat/mounting/pi-henkes/self-alignment/conf_PBC/long/","PBC_0*.dat")

fp_ra = scan_dat("/Users/kammeraat/mounting/pi-henkes/self-alignment/conf_PBC/long/","last*.dat")



frames_name, frames_data = read_frames(fp_sim)

support_name, frames_support = read_frames(fp_ra)

Nint = length(frames_data[1].x)
Nt =  length(frames_data)

x = zeros(Nint, Nt)
y = zeros(Nint, Nt)

vx = zeros(Nint, Nt)
vy = zeros(Nint, Nt)

R = collect(frames_data[1].radius)
type = collect(frames_data[1].type)
@views for i in 1:Nt
    x[:,i] .=  frames_data[i].x
    y[:,i] .=  frames_data[i].y
    
    vx[:,i] .=  frames_data[i].vx
    vy[:,i] .=  frames_data[i].vy
end

x0 = collect(frames_support[end].x)

y0 = collect(frames_support[end].y)


include("../analysis/AnalysisPipeline.jl") 


D = construct_D(x0, y0, 1, R, type, periodic_system_sizes = [50., 50., 0.])


eigenmodes = diagonalize_D(D)

v_projs = project_on_eigvecs(eigenmodes["eigvecs"], vx , vy)

eigval_pos  =eigenmodes["eigvals"].>0

the_eigvals = eigenmodes["eigvals"][eigval_pos]

tau = 20

theory_ABP = @. 1 ./(2 .+ 2*tau .* the_eigvals)


begin
v0 = 0.01
Dr = 1/20
J=0
v_projs2 = mean(v_projs.^2, dims=2)



f = Figure()
ax = Axis(f[1,1], title="Dr = $Dr, J=$J", xlabel="λ", ylabel=L"|\dot{a}_ν(t)|^2/v_0^2", yscale=log10, xscale=log10)



scatter!(ax,the_eigvals, v_projs2[eigval_pos])# ,  label="J = $(e["J"])",alpha=0.3)

        #lines!(ax,the_eigvals,theory_ABP, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"]) ABP theory ", alpha=0.2)
       
        #scatterlines!(ax,the_eigvals,real.(B), color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"])", linestyle=:dash)

        
#lines!(ax,the_eigvals,theory_integral, linestyle=:solid)#,  label="J = $(e["J"])")
lines!(ax,the_eigvals,theory_ABP.*v0^2, linestyle=:solid, color=:red)#

xlims!(1e-3,1e1)
display(f)

end



base_folder  = "/Users/kammeraat/compare_SAMoSA_ABP/"
rp, sp, ap = auto_analysis_dir(base_folder, "sa_raw_data.h5"; support_raw_data_file_name_pattern = "ra_raw_data.h5")

run_sequential_analysis(rp, ap,run_free_BP_analysis!, support_raw_data_file_paths=sp, overwrite=true)

JAMs = load_file("/Users/kammeraat/compare_SAMoSA_ABP/phi_1.3/N_800/analysis/v0_0.01.h5")
using CairoMakie
begin
v0 = 0.01
Dr = 1/20
J=0
v_projs2 = mean(v_projs.^2, dims=2)


JAMS_v_projs2 = mean( JAMs["projs"]["v_projs"].^2, dims=2)

JAMS_v_projs2_eigvals = JAMs["eigenmodes"]["eigvals"]

JAMS_v_projs2_eigvals_pos = JAMS_v_projs2_eigvals .>0


f = Figure()
ax = Axis(f[1,1], title="Dr = $Dr, J=$J", xlabel=L"λ_{\nu}", ylabel=L"|\dot{a}_ν(t)|^2", yscale=log10, xscale=log10)



scatter!(ax,the_eigvals, v_projs2[eigval_pos], label="SAMoS")# ,  label="J = $(e["J"])",alpha=0.3)

        #lines!(ax,the_eigvals,theory_ABP, color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"]) ABP theory ", alpha=0.2)
       
        #scatterlines!(ax,the_eigvals,real.(B), color=e["J"], colorrange = (0, .1) ,  label="J = $(e["J"])", linestyle=:dash)

scatter!(ax,JAMS_v_projs2_eigvals[JAMS_v_projs2_eigvals_pos], JAMS_v_projs2[JAMS_v_projs2_eigvals_pos],label="JAMs")
#lines!(ax,the_eigvals,theory_integral, linestyle=:solid)#,  label="J = $(e["J"])")
lines!(ax,the_eigvals,theory_ABP.*v0^2, linestyle=:solid, color=:red)#
f[1,2]=Legend(f,ax)
xlims!(1e-3,1e1)
display(f)
save("/Users/kammeraat/JAMs_SAMoS_vel_proj2.pdf",f,backend=CairoMakie)

end