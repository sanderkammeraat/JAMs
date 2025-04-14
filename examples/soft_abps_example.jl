include(joinpath("..","src","Engine.jl"))
function simulation()

    external_forces = [ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1)]

    pair_forces = [soft_disk_force(1,2)]

    
    
    local_dofevolvers =  [overdamped_θω_evolver(1), overdamped_xvf_evolver(1)]
    global_dofevolvers = []

    #Initialize state
    N=100
    ϕ = 0.5
    L=sqrt(N*pi/ϕ)
    poly = 15e-2
    seed = 1
    Random.seed!(seed)
    initial_state = [ PolarParticle2d([i],[1],[1],[0.3],[0.01],rand(Uniform(-L/2, L/2),2),[0,0],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],rand(Uniform(1-poly, 1+poly),1),[1.],[0.,0.],[0.,0.],[0,0]) for i=1:N];

    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, true,10.)#2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.1,1e4,Tplot=10, seed=2, Tsave = nothing, save_folder_path=joinpath("/Users/kammeraat/mounting/data1_kammeraat/testing/network", "prof_soft_abps"), save_functions=[save_2d_polar_θ!], fps=120, plot_functions=[plot_disks!, plot_directors!, plot_velocity_vectors!]);
    return sim
end
#
sim = simulation(); 

@profview simulation()

@time simulation();


file = jldopen(joinpath(pwd(), "prof_soft_abps","raw_data.jld2"),"r")

file["integration_info"]["master_seed"]
file["frames"]
file
close(file)

struct te
    r::Float64

end

test = te(1)

hasproperty(test, :s)