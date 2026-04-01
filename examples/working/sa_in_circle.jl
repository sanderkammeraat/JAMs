
include(joinpath("../../src","Engine.jl"))

function simulation()



    pair_forces = (soft_disk_force([1,2],[1 .05; .05 1]),)
    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_force(1,100),ABP_perpendicular_angular_noise(1,[0,0,1]))

    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1),)
    global_dofevolvers = ()
    field_dofevolvers = []
    N = 2
    r_avg = 0.3
    R = 2.2
    bdens = 2.
    Nb = floor(Int64,bdens*2*pi*R/(2*r_avg))

    r_max = r_avg 

    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[ PolarParticle3d([i], [1], [1.], [1.], [1], [0.2], [0.01], [rand(Uniform(-2*R/3, 2*R/3)) , rand(Uniform(-2*R/3,2*R/3)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];
  
    j=1
    for i in N+1:N+Nb

        push!(initial_state,ConfinedPolarParticle3d([i], [2], [1.],[1.], [0.2], [0.], [0.0], [R*cos(2*pi/Nb*j) , R*sin(2*pi/Nb*j),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],[0,0,1],[0,0,0],[0,0,0]))
        j+=1
    end
    Lx = 2*R+2*r_max
    Ly = 2*R+2*r_max

    sizes = [Lx,Ly,10];
    print(sizes)
    initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,2.5*2);

    #Run integration
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting
    sim = Euler_integrator(system,0.001, 1e4, Tplot=1,fps=60,plot_functions=(plot_disks_orientation!,plot_directors!), plotdim=2, Tsave=nothing, save_folder_path = "/Users/kammeraat/test_saving_speed/", save_functions=[save_2d_polar_p!]); 
    return sim;

end

simulation()