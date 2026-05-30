include(joinpath("..","..","src","Engine.jl"))

function simulation()

    
    external_forces = (Lorenz_system(ontypes=1,σ = 10, β = 8/3, ρ = 28),)
    pair_forces = []
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),)
    global_dofevolvers = []
    field_dofevolvers = []

    N=1000
    ϕ = 0.5
    poly=15e-2
    Rs = rand(Uniform(1-poly, 1+poly),N)
    display(size(Rs))

    L =  150
    Ly = L
    Lx = L
    Lz = L
    #rand(Uniform(-1,1),3)
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [0.3], [0.001],[1,1,1] .+5*rand(Uniform(-1,1),3) ,[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];


    display(L)
    sizes = [Lx,Ly,Lz];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, false,2.5*(1+poly));

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.0001, 1e4,sbs=true, Tplot=100,fps=60,plot_functions=(plot_trajectories!,plot_sized_points!), plotdim=3, Tsave=nothing); 
    return sim;

end


sim = simulation()  