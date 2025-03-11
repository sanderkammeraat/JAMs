include(joinpath(pwd(),"src","Engine.jl"))
function simulation(seed,save_folder_path)

    external_forces = [ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1)]

    pair_forces = [soft_disk_force(1,1)]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1000
    ϕ = 0.9
    L=sqrt(N*pi/ϕ)
    poly = 15e-2
    init_seed = 1
    Random.seed!(init_seed)
    initial_state = [ PolarParticle2d(i,1,1,0.3,0.1,rand(Uniform(-L/2, L/2),2),[0,0],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.],[0,0]) for i=1:N];

    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.01,20,seed=seed,Tsave=100,save_functions=[save_2d_polar_θ!], save_folder_path=save_folder_path)#, Tplot = 100, fps=120, plot_functions=[plot_disks!]);
    return sim
end
#

begin
seeds = [1, 2, nothing, 1, 10000, 10000,1]
save_folder_paths = []

for (i, seed) in pairs(seeds)

    save_folder_path_i = joinpath(pwd(),"testing","setting_seeds","t_"*string(Threads.nthreads())*"_i_"*string(i))

    simulation(seed,save_folder_path_i); 

    push!(save_folder_paths,save_folder_path_i)

end

xs = zeros(length(seeds))
for (i, save_folder_path) in pairs(save_folder_paths)

    file = jldopen(joinpath(save_folder_path,"raw_data.jld2"), "r")
    xs[i] = file["frames"]["10"]["x"][4]
end


@assert xs[1]!=xs[2]

@assert xs[1]!=xs[3]

@assert xs[1]==xs[4]==xs[7]


end

#Rerun previous code block with 4 threads and then run the 

x1s = zeros(length(seeds))
x4s = zeros(length(seeds))
for (i, save_folder_path) in pairs(seeds)

    file_t1 = jldopen(joinpath(pwd(),"testing","setting_seeds","t_"*"1"*"_i_"*string(i),"raw_data.jld2"), "r")
    file_t4 = jldopen(joinpath(pwd(),"testing","setting_seeds","t_"*"4"*"_i_"*string(i),"raw_data.jld2"), "r")

    x1s[i] = file_t1["frames"]["10"]["x"][4]

    x4s[i] = file_t4["frames"]["10"]["x"][4]
end

@assert x1s[1]==x4s[1]

@assert x1s[2]==x4s[2]

@assert x1s[3]!=x4s[2]

print(x1s)
print(x4s)

#Test seed=nothing reproducibility

f3 = jldopen(joinpath(pwd(),"testing","setting_seeds","t_"*"1"*"_i_"*string(3),"raw_data.jld2"), "r")
seed_rep = f3["integration_info"]["master_seed"]

seed_rep_jl = jldopen(joinpath(pwd(),"testing","setting_seeds","t_"*"1"*"_i_"*string(3),"JAMs_container.jld2"), "r")["integration_info"]["master_seed"]

@assert seed_rep == seed_rep_jl

save_folder_path_3 = joinpath(pwd(),"testing","setting_seeds","c_"*"1"*"_i_"*string(3))

simulation(seed_rep,save_folder_path_3); 

c3 = jldopen(joinpath(save_folder_path_3,"raw_data.jld2"),"r")


f3["frames"]["10"]["x"][4]
c3["frames"]["10"]["x"][4]
@assert f3["frames"]["10"]["x"][4] == c3["frames"]["10"]["x"][4]



