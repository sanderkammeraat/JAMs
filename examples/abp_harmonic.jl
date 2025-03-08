include(joinpath("..","src","Engine.jl"))

function simulation()

    external_forces = (ABP_3d_propulsion_force(1), ABP_3d_angular_noise(1),external_harmonic_force(1,0.1),self_align_with_v_force(1,0.1))

    pair_forces = []#[soft_disk_force(1,1)]

    
    
    dofevolvers =  [overdamped_evolver!]

    #Initialize state
    N=1
    L=40.

    poly = 1e-6

    
    initial_state = [ PolarParticle3d(i,1,1.,1.,rand(Uniform(1-poly, 1+poly)),2.,0.001, rand(Uniform(-L/2, L/2),3),[0.,0.,0.],[0.,0.,0.],[0.,0.,0.],[0.,0.,0.],normalize(rand(Normal(0, 1),3)),[0,0,0],[0,0,0]) for i=1:N];
    
    size = [L,L,L];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,4.);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system, 0.01, 1000, Tplot=10, fps=120,plot_functions=(plot_sized_points!, plot_directors!));
    return sim

end


sim = simulation();

