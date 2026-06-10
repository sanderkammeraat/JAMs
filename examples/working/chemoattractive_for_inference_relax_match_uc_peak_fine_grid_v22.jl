

include(joinpath("..","..","src","Engine.jl"))

#sim units
# [t] = s
# [L] =  um

function simulation(save_path)

    external_forces = (ABP_perpendicular_angular_noise(1,[0,0,1]),)
    R = 5.
    pair_forces =(exp_repulsion_force(1, 5. /R*R , 0.45),)
    #type,  consumption, cmid, v0max, σ in Lorentzian
    #asymmetric_field_propulsion_distr_force(1,2000,log10(0.2),4.5*20,0.261)
    #c0 = 0.0005
    #cmid = 0.001

    

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers =  (overdamped_xvf_evolver(1), overdamped_pq_xyc_evolver(1))
    global_dofevolvers =  ()
    field_dofevolvers = (overdamped_CCvCf_evolver(1),)

    ϕ = 0.01
    L =  300. *R #300

    N = floor(Int64, L^2 * ϕ / ( π * R^2))
    display(N)
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [R], [0.3], [0.015], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];

    Lx=L
    Ly=L
    sizes = (Lx,Ly,2.);
    lbin = 1.
    print(lbin)

    z_bin_centers = [0.]
    #Consumption is tied to dirac delta, which has inverse bin size^2 dimensions which need to be compensated
    # One R comes from lbin, the other for scaling the magnitudes of the active forces due to the speed, which should not affect the torque
    field_forces =(field_propulsion_distr_force(1,5. *R^2 / lbin^2,log10(2),4.5*2.5,0.261), self_align_with_∇C_force(1,3 * R), grad_field_propulsion_force(1,0.,-10 *R^2))
    
    x_bin_centers = [-Lx-lbin]
    x_bin_centers = append!(x_bin_centers,range(start=-Lx/2, stop=Lx/2+0.1*lbin, step=lbin).+lbin/2)
    x_bin_centers = append!(x_bin_centers,Lx+lbin)

    y_bin_centers = [-Ly-lbin]
    y_bin_centers = append!(y_bin_centers,range(start=-Ly/2, stop=Ly/2+0.1*lbin, step=lbin).+lbin/2)
    y_bin_centers = append!(y_bin_centers,Ly+lbin)
    bin_centers = [x_bin_centers, y_bin_centers,z_bin_centers]

    C = 1 .*ones(length(x_bin_centers), length(y_bin_centers))
    
    initial_field_state=[GeneralField2d(1,1,bin_centers,C, C.*0, C.*0, lbin)]
    # [field_propulsion_3d_force(1,1e-2,0.01)]
    field_updaters = (Diffusion(1,R^2*1e1),Relax(1,1.,0.3), GhostSet(1))


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers, global_dofevolvers, field_dofevolvers, true,5*R);

    #Run integration
    #Use plot_disks! for nice visuals
    #Use plot_points! for fast plotting
    #Tsave 20

    sim = Euler_integrator(system,0.0005 ,600, Tplot=nothing, fps=120, res=(1000,1000), Tsave = nothing ,save_folder_path = save_path, save_functions=(save_2d_polar_p!, save_single_2d_field_10!),plot_functions=(plot_field_magnitude!, plot_disks_orientation!, plot_directors!), plotdim=2); 
    return sim

end

sim = simulation(joinpath(homedir(),"surfdrive","ActivePolygonClusters","simulations","for_inference_v23_fine_grid","phi_0p01","simdata"))


#sim = simulation(joinpath("/Volumes","T7_Shield","ActivePolygonClusters","simulations","for_inference_v23_fine_grid","phi_0p01","simdata"))


@profview  simulation(joinpath(homedir(),"surfdrive","ActivePolygonClusters","simulations","for_inference_v23_exp_rep","phi_0p01","simdata"))



