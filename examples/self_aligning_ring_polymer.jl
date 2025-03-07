include("../src/Engine.jl")

function simulation()

    external_forces = [ ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,0.4),ABP_perpendicular_angular_noise(1,[0,0,1])]

    

    #dofevolvers = [inertial_evolver!]
    dofevolvers = [overdamped_evolver!]
    N=100
    
    l=2.
    k=5.
    R=N*l/2/pi
    L=5*R
    v0=0.5
    Dr=0.01
    pair_forces = [periodic_chain_force(1,k,l,1,N),soft_disk_force([1],1e0)]
    initial_state = [PolarParticle3d(i, 1,1, 1, 1, v0, Dr, [R*cos(2*pi/N*i) , R*sin(2*pi/N*i),0],[0,0,0],[0.,0.,0.], [0,0,0],[0,0,0],push!(normalize(rand(Normal(0),2)),0).*[1,1,0],[0,0,0],[0,0,0]) for i in 1:N];

    #push!(initial_state, NewPolarParticle3d(N, 1, 1, 1, 0, 0.01, [0 , -N*2,0],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0]))
    
    size = [L,L,10.];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,10.);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,1e-2, 1e5, Tplot=5e1, fps=120, plot_functions=(plot_sized_points!, plot_directors!, plot_velocity_vectors!), plotdim= 2); 
    return sim

end


sim = simulation()

