include("../src/Engine.jl")
using Random, Distributions

function simulation()

    external_forces = [ABP_2d_angular_noise(1)]

    pair_forces = [soft_disk_force(1,1)]

    field_forces = [field_propulsion_force(1,1e-2,0.01)]
    field_updaters = [PeriodicDiffusion(1,8e-3)]

    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1000
    Lx=120.
    Ly=120.
    L=min(Lx,Ly)
    poly = 0.0000002


    initial_particle_state = [ PolarParticle2d(i,1,1,0.0,0.0001,[rand(Uniform(-L/2, L/2)) ,rand(Uniform(-L/2,L/2))],[0.,0.],[0.,0.],[rand(Uniform(-pi, pi))],[0.],1,rand(Uniform(1-poly, 1+poly)),[0.,0.],[0.,0.]) for i=1:N];
    
    size = [Lx,Ly];

    lbin = 1
    print(lbin)

    x_bin_centers = range(start=-Lx/2, stop=Lx/2, step=lbin).+lbin/2
    y_bin_centers = range(start=-Ly/2, stop=Ly/2, step=lbin).+lbin/2
    bin_centers = [x_bin_centers, y_bin_centers]
    print(x_bin_centers)

    v00=0.4

    C = ones(length(x_bin_centers), length(y_bin_centers))*v00
    
    initial_field_state=[FuelField2d(1,bin_centers,C, C.*0, C.*0)]

    system = System(size, initial_particle_state,initial_field_state,external_forces, pair_forces, field_forces, field_updaters,dofevolvers, true, 1e1);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    plot_functions = (plot_sized_points!, plot_directors!, plot_velocity_vectors!,plot_field_magnitude!)
    sim = Euler_integrator(system, 0.1, 1e7, 1e7,5, 120, plot_functions,2);
    #particle_states,field_states = Euler_integrator(system, 0.1, 10, 10000000000,0, 0, plot_functions,false);
    return sim

end


sim = simulation()

@time simulation();





psO = Observable(particle_states[1])



f, ax = setup_system_plotting(system.sizes,[test!], false, psO,[],0)

ax=test!(ax, psO ,[])

fps = 1
for i in eachindex(particle_states)
    psO[].= particle_states[i]
    notify(psO)
    sleep(1/fps)
end

GLMakie.activate!()
f = Figure()