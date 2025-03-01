include("../src/Engine.jl")

#%%Running
using Random, Distributions

function simulation()

    external_forces = (ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1))

    pair_forces = []#[soft_disk_force(1,1)]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1000
    ϕ = 0.2
    L=sqrt(N*pi/ϕ)
    poly = 0.2


    
    initial_state = []  
    for i=1:N
    xi = [rand(Uniform(-L/2, L/2)) ,rand(Uniform(-L/2, L/2))]
    push!(initial_state,PolarParticle2dSave(i,1,1,0.3,0.1,xi,xi,[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.],[0,0]));
    end

    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,4e0);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.1, 1000, 100, 5,120);
    return sim
end


sim= simulation();

make_movie(sim,"/Users/kammeraat/test_JAMS/movies/soft_abps_v21.mp4",(plot_disks!,plot_directors!, plot_velocity_vectors!),60)



folder_path = "/Users/kammeraat/test_JAMS/savedata/"

file_path=save_SIM(folder_path, "abps_sim_Dr", sim)


siml = load_SIM(file_path);

