include(joinpath("..","..","src","Engine.jl"))


function soft_disk_no_overlap()
    pair_forces =[soft_disk_force([1,2],[1 1; 1 1])]#[soft_shape_disk_force(1,0)]

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver([1,2]),overdamped_2d_shape_evolver([1,2]))
    global_dofevolvers = []
    field_dofevolvers = []
    N1=150
    N2 = 150
    N= N1+N2
    ϕ =0.5
    Rno=3
    L =  sqrt(pi * N*Rno^2 / ϕ)
    xo = Float64[ -1 0 0;  1/2 -1/2*sqrt(3)   0 ;  1/2 1/2*sqrt(3)   0]

    initial_state = PolarShape[PolarShape([i],[1], [1], [1], [Rno], [0.3], [0.001], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0],deepcopy(xo),xo,[1,1,1]) for i=1:N1 ]


    xo_inv = Float64[ 1 0 0;  -1/2 -1/2*sqrt(3)   0 ;  -1/2 1/2*sqrt(3)   0]
    for i=N1+1:N2+N1
        push!(initial_state, PolarShape([i],[2], [1], [1], [Rno], [0.3], [0.001], [rand(Uniform(-L/2, L/2)) , rand(Uniform(-L/2,L/2)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0],deepcopy(xo_inv),xo_inv,1*ones(size(xo_inv)[1])))
    end
    sizes = [L,L,4];
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    external_forces = []

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true, 10.);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.1, 50, Tplot=nothing,fps=120,plot_functions=(plot_shape_disks!,plot_directors!), plotdim=2); 
    return sim;

end
function simulation(soft_disk_no_overlap_result)

    

    pair_forces =[exp_shape_disk_force([1,2],[1 1 ; 1 1])]

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver([1,2]),overdamped_pq_evolver([1,2]),overdamped_2d_shape_evolver([1,2]))
    global_dofevolvers = []
    field_dofevolvers = [overdamped_CCvCf_evolver(1)]

    initial_particle_state =deepcopy(soft_disk_no_overlap_result.final_particle_state)

    #Modify initial state
    for (i, p_i) in pairs(initial_particle_state)
        #reinitialize forces 
        p_i.f.*=0
        p_i.q.*=0
        p_i.R.=1

    end
    sizes = soft_disk_no_overlap_result.system.sizes;

    Lx=sizes[1]
    Ly=sizes[2]
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
    
    initial_field_state=[FuelField2d(1,1,bin_centers,C, C.*0, C.*0)]
    # [field_propulsion_3d_force(1,1e-2,0.01)]
    field_updaters = (PeriodicDiffusion(1,1e-1),AvgSetwoGhost(1,1.), GhostSet(1))
    field_forces = (field_propulsion_distr_force(1,0.2,1.,0.5,1),self_align_with_∇C_force(1,-1))

    external_forces = (ABP_perpendicular_angular_noise([1,2],[0,0,1]), )

    system = System(sizes, initial_particle_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true, 10.);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.01, 1e4, Tplot=40,fps=120,plot_functions=(plot_shape_disks_type!,plot_disks_orientation!,plot_directors!, plot_field_magnitude!), plotdim=2); 
    return sim;

end

soft_disk_no_overlap_result = soft_disk_no_overlap()
sim = simulation(soft_disk_no_overlap_result)  

    