include(joinpath("..","src","Engine.jl"))

#Assuming that this folder contains an abp simulation
folder_path = joinpath(pwd(), "testing", "correctness_rng")
file_path = joinpath(folder_path,"raw_data.jld2")

f = jldopen(file_path, "r")


#collect positions of particles
Nt = length(f["frames"])
Np = length(f["frames"]["1"]["xuw"])
r = zeros(Nt,Np, 2);

@views for frame in eachindex(f["frames"])

    r[parse(Int, frame),:,1] = f["frames"][frame]["xuw"]

    r[parse(Int, frame),:,2] = f["frames"][frame]["yuw"]

end

size(r)[1]
@views function calculate_MSD(r)

    Nt = size(r)[1]

    msd = zeros(Nt)

    for i in 1:Nt

        

        msd[i] = mean( (r[i:end,:,1]-r[1:end-i+1,:,1]) .^2 + (r[i:end,:,2]-r[1:end-i+1,:,2]) .^2 )


    end
    return msd
end

msd = calculate_MSD(r)
msd[3]  

v0 = f["frames"]["1"]["v0"][1]
Dr = f["frames"]["1"]["Dr"][1]

begin
fig= Figure()
ax = Axis(fig[1,1],xscale=log10,yscale=log10, xlabel="Δt", ylabel="msd")
t = f["integration_info"]["save_tax"]

numplot=plot!(ax,t,msd)

theplot=lines!(t, 2*v0^2/Dr *(t -1/Dr*(1 .-exp.(-t*Dr))), color=:tomato)

slope1line = lines!(t, 200*t.*2*v0^2*Dr, color=:black, linestyle=:dash)

slope2line = lines!(t[t.<1e1], t[t.<1e1].^2, color=:black, linestyle=:dot)

xlims!(1e-2,1e4)

ylims!(1e-2,1e4)

Legend(fig[1,2],[numplot, theplot, slope1line,slope2line], ["JAMs numerical", "2*v0^2/Dr *(t -1/Dr*(1 .-exp.(-t*Dr)))", "Slope 1", "Slope 2"])

display(fig)
save(joinpath(folder_path,"msd.png"),fig)
end 

