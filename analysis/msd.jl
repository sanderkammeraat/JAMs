include("../src/Engine.jl")


file_path = "/Users/kammeraat/test_JAMS/savedata/abps_sim_large_long.jld2"

siml = load_SIM(file_path);

particle_states = siml.particle_states

#collect positions of particles
Nt = length(particle_states)
Np = length(particle_states[1])
r = zeros(Nt,Np, 2)

@views for i in eachindex(particle_states)

    for j in 1:Np
        r[i,j,:].=particle_states[i][j].xs


    end
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


f= Figure()
ax = Axis(f[1,1],xscale=log10,yscale=log10)
ylims!(ax,1e0,1e3)

plot!(ax,range(0,length(msd)-1),msd)

display(f)
