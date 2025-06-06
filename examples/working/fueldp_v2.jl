include(joinpath("..","..","src","Engine.jl"))

function simulation()

    external_forces = [ABP_perpendicular_angular_noise(1,[0,0,1])]

    pair_forces = [soft_disk_force(1,1.)]

    field_forces =[asymmetric_field_propulsion_distr_force(1,0.001,1.,0.3,10), self_align_with_∇C_unit_force(1,-0.01)]

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers =  [overdamped_xvf_evolver(1), overdamped_pq_evolver(1)]
    global_dofevolvers =  []
    field_dofevolvers = [overdamped_CCvCf_evolver(1)]
    N=1000
    ϕ = 0.1
    L =  sqrt(N *  π * 1^2 / ϕ)
    compres = 1.
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [1], [0.3], [0.01], [rand(Uniform(-L/2*compres, L/2*compres)) , rand(Uniform(-L*compres/2,L*compres/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    Lx=L
    Ly=L
    size = [Lx,Ly,2];
    lbin = sqrt(2)
    print(lbin)

    z_bin_centers = [0.]
    x_bin_centers = [-Lx-lbin]
    x_bin_centers = append!(x_bin_centers,range(start=-Lx/2, stop=Lx/2, step=lbin).+lbin/2)
    x_bin_centers = append!(x_bin_centers,Lx+lbin)

    y_bin_centers = [-Ly-lbin]
    y_bin_centers = append!(y_bin_centers,range(start=-Ly/2, stop=Ly/2, step=lbin).+lbin/2)
    y_bin_centers = append!(y_bin_centers,Ly+lbin)
    bin_centers = [x_bin_centers, y_bin_centers,z_bin_centers]

    C = ones(length(x_bin_centers), length(y_bin_centers))
    
    initial_field_state=[FuelField2d(1,1,bin_centers,C, C.*0, C.*0)]
    # [field_propulsion_3d_force(1,1e-2,0.01)]
    field_updaters = [PeriodicDiffusion(1,1), GhostSet(1)]




    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, true,2.5);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,5e-2, 6e2,Tsave=nothing, save_folder_path=joinpath(homedir(), "ADCA","simulations","anti-align","simdata"), save_functions=(save_2d_polar_p!,save_single_2d_field!), Tplot= 2e1, fps=120, plot_functions=(plot_field_magnitude!,plot_disks!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end


sim = simulation()