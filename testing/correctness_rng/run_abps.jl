include(joinpath(pwd(),"src","Engine.jl"))
function simulation(save_folder_path)

    external_forces = [ABP_2d_propulsion_force(1), ABP_2d_angular_noise(1)]

    pair_forces = []

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1000
    ϕ = 0.9
    L=sqrt(N*pi/ϕ)
    poly = 15e-2
    initial_state = [ PolarParticle2d(i,1,1,0.3,0.1,rand(Uniform(-L/2, L/2),2),[0,0],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],rand(Uniform(1-poly, 1+poly)),1.,[0.,0.],[0.,0.],[0,0]) for i=1:N];

    size = [L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.01,100,Tsave=100,save_functions=[save_2d_polar_θ!], save_folder_path=save_folder_path)#, Tplot = 100, fps=120, plot_functions=[plot_disks!]);
    return sim
end
#

save_folder_path = joinpath(pwd(),"testing","correctness_rng")

sim = simulation(save_folder_path); 
