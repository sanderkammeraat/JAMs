include(joinpath("..","..","src","Engine.jl"))

function simulation(save_path)

    external_forces = [ABP_perpendicular_angular_noise(1,[0,0,1])]

    pair_forces = [soft_disk_force(1,1.)]
    #type,  consumption, cmid, v0max, σ in Lorentzian
    field_forces =(field_propulsion_distr_force(1,0.2, log10(1) ,10 ,30), self_align_with_∇C_force(1,5), grad_field_propulsion_force(1,0.,-2))
    #dofevolvers = [inertial_evolver!]
    local_dofevolvers =  [overdamped_xvf_evolver(1), overdamped_pq_xyc_evolver(1)]
    global_dofevolvers =  []
    field_dofevolvers = [overdamped_CCvCf_evolver(1)]
    N=1
    ϕ = 0.1
    L =  10.
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [1], [0.3], [0.001], [0.,L/3,0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    Lx=L*3
    Ly=L*3
    sizes = @SVector [Lx,Ly,2];
    lbin = 1
    print(lbin)

    z_bin_centers = [0.]
    x_bin_centers = [-Lx-lbin]
    x_bin_centers = append!(x_bin_centers,range(start=-Lx/2, stop=Lx/2+0.1*lbin, step=lbin).+lbin/2)
    x_bin_centers = append!(x_bin_centers,Lx+lbin)

    y_bin_centers = [-Ly-lbin]
    y_bin_centers = append!(y_bin_centers,range(start=-Ly/2, stop=Ly/2+0.1*lbin, step=lbin).+lbin/2)
    y_bin_centers = append!(y_bin_centers,Ly+lbin)
    bin_centers = [x_bin_centers, y_bin_centers,z_bin_centers]

    C = ones(length(x_bin_centers), length(y_bin_centers))

    midinds = [round(Int64,length(x_bin_centers)/2), round(Int64, length(y_bin_centers)/2)]

    display(midinds)
    display(size(C))
    C[midinds...] = 0.




    #C = (C .+ reverse(C, dims=1)) ./2
    
    initial_field_state=[FuelField2d(1,1,bin_centers,C, C.*0, C.*0)]
    # [field_propulsion_3d_force(1,1e-2,0.01)]
    field_updaters = [PeriodicDiffusion(1,1e0), AvgSetwoGhost(1,1.),GhostSet(1), IndSet(1, midinds, 0.)]


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, true,2.5);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.05, 1e4, Tplot= 20, fps=120, Tsave = nothing ,save_folder_path = save_path, save_functions=(save_2d_polar_p!,save_single_2d_field!),plot_functions=(plot_field_magnitude!,plot_disks_orientation!, plot_directors!), plotdim=2); 
    return sim

end


sim = simulation(joinpath(homedir(),"surfdrive","ActivePolygonClusters","simulations","patch_experiment","simdata"))
