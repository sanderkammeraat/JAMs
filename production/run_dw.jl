include(joinpath("..","src","Engine.jl"))

function put_in_well_a_or_b(i,N,xa,xb)

    xi = i < N/2 ? xa : xb
    return xi

end

function simulation(k, Dr, J, Tplot, Tsave,tend)

    

    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    pair_forces =()#(soft_disk_force(1,0),)


    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = ()
    field_dofevolvers = ()

    #Number of particles, 1000 in each well
    N=2000

    Lx =  10.

    Ly = Lx

    

    xa = -1
    xb = 1

    ka = k
    kb = ka

    sigma = .2
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [1.], [0.3], [Dr], [put_in_well_a_or_b(i,N,xa,xb) + rand(Uniform(-sigma,sigma)) , rand(Uniform(-sigma,sigma)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];



    sizes = (Lx,Ly,4.);
    initial_field_state=[]
    field_forces = ()
    field_updaters = ()

    #β=-1 interesting!
    external_forces = (external_double_gaussian_force(1,ka,kb, [xa,0,0], [xb,0,0]),ABP_perpendicular_angular_noise(1,[0,0,1]),self_align_with_v_unit_force(1,J),ABP_3d_propulsion_force(1))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, false,3.);

    #Run integrationov
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting

    #save_folder = "/Users/kammeraat/dwsa/single/simdata/v0_$v0/Dr_$Dr/J_$J/"
    save_folder = "/Volumes/T7_Shield/sa_double_well/vary_k/k_$k/"
    sim = Euler_integrator(system,0.01,tend, Tsave=Tsave, fps=120,Tplot=Tplot,plot_functions=(plot_potential!,plot_trajectories!), plotdim=2, save_folder_path = save_folder, save_functions = (save_2d_polar_p!,),res=(1000,1000)); 
    return sim;

end

#v0, Dr, J, Tplot, Tsave
#start with k=0.1, k=2 is max, close to 1.55 we get hopping, regular 2- orbits round 0.3

for k in [0.1, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.5, 1.55, 1.6, 1.8, 2.0]
    #k =1.6
    display(k)
    Dr = 0.01
    J = 1
    Tplot =nothing #plot every nth timestep, set to Tplot=nothing to turn off plotting

    Tsave = 10 #save every nth timestep, set to Tsave=nothing to turn off saving

    tend = 1e3#1e3 #for how long to run the simulation (in units of simulation time)

    simulation(k, Dr ,J, Tplot, Tsave,tend)
end