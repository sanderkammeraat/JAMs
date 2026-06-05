include(joinpath("..","..","src","Engine.jl"))

#sim units
# [t] = s
# [L] = 5 um

# Then v0 = 0.4 corresponds to 0.4 * 5 = 2 um/s
#.4, 0.8
function simulation(save_path)

    external_forces = (ABP_perpendicular_angular_noise(1,[0,0,1]),)
    R = 1.
    pair_forces =(exp_repulsion_force(1, 5. *R,0.5),)
    #type,  consumption, cmid, v0max, σ in Lorentzian
    #asymmetric_field_propulsion_distr_force(1,2000,log10(0.2),4.5*20,0.261)
    
    field_forces =(asymmetric_field_propulsion_distr_force(1,1000,1.,1.5/0.3*5,0.3), self_align_with_∇C_force(1,5), grad_field_propulsion_force(1,0.,-15))

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers =  (overdamped_xvf_evolver(1), overdamped_pq_xyc_evolver(1))
    global_dofevolvers =  []
    field_dofevolvers = (overdamped_CCvCf_evolver(1),)

    ϕ = 0.1
    L =  300*R #300

    N = floor(Int64, L^2 * ϕ / ( π * R^2))
    display(N)
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [R], [0.3], [0.015], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    Lx=L
    Ly=L
    size = [Lx,Ly,2];
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

    C = 1 .*ones(length(x_bin_centers), length(y_bin_centers))
    
    initial_field_state=[FuelField2d(1,1,bin_centers,C, C.*0, C.*0)]
    # [field_propulsion_3d_force(1,1e-2,0.01)]
    field_updaters = (Diffusion(1,1e1),Relax(1,1.,0.3), GhostSet(1))


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, true,5*R);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    #Tsave 20

    sim = Euler_integrator(system,0.005, 1e4, Tplot= 20, fps=120, res=(1000,1000), Tsave = nothing ,save_folder_path = save_path, save_functions=(save_2d_polar_p!,save_single_2d_field!),plot_functions=(plot_field_magnitude!, plot_disks_orientation!, plot_directors!), plotdim=2); 
    return sim

end

sim = simulation(joinpath(homedir(),"surfdrive","ActivePolygonClusters","simulations","for_inference_v17_exp_rep","phi_0p01","simdata"))

@profview_allocs simulation(joinpath(homedir(),"surfdrive","ActivePolygonClusters","simulations","for_inference_v17_exp_rep","phi_0p01","simdata"))
