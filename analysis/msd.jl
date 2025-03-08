include(joinpath("..","src","Engine.jl"))


file_path = "/Users/kammeraat/test_JAMS/msd_analysis/raw_data.jld2"

f = jldopen(file_path, "r")


#collect positions of particles
Nt = length(f["frames"])
Np = length(f["frames"]["1"]["xuw"])
r = zeros(Nt,Np, 2)

@views for frame in eachindex(f["frames"])

    r[parse(Int, frame),:,1] = f["frames"][frame]["xuw"]

    r[parse(Int, frame),:,2] = f["frames"][frame]["yuw"]

end

size(r)[1]
@views function calculate_MSD(r)

    Nt = size(r)[1]

    msd = zeros(Nt)

    for i in 1:Nt

        msd[i] = mean((r[i+1:end,:,:].-r[1:end-i,:,:]).^2)


    end
    return msd
end

msd = calculate_MSD(r)
msd[3]

begin
fig= Figure()
ax = Axis(fig[1,1],xscale=log10,yscale=log10)

plot!(ax,range(0,length(msd)-1)*f["integration_info"]["Tsave"]*f["integration_info"]["dt"],msd)

display(fig)
end