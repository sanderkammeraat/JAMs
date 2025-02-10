
include("../src/Engine.jl")

#%%Running
using Random, Distributions

function simulation()

    external_forces = (ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1))

    pair_forces = [soft_disk_force(1,1)]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=500
    ϕ = 0.2
    L=sqrt(N*pi/ϕ)
    poly = 0.2


    initial_state = [ PolarParticle2d(i,1,1,0.3,0.0001,[rand(Uniform(-L/2, L/2)) ,rand(Uniform(-L/2,L/2))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.]) for i=1:N];


    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,4e0);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.1, 20, 10, 5,120);
    return sim
end


sim= simulation();


folder_path = " folder path ending in  /"

file_path=save_SIM(folder_path, "soft_abps_sim", sim)


siml = load_SIM(file_path);

