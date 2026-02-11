include(joinpath("..","..","src","Engine.jl"))


function simulation()

    L = 10   
    

    phi = 1.0      #packing fraction; 0.005 for 1 particle
    dt = 0.0001              #timestep
    tot_time = 10
    #steps = int(tot_time/dt)
                 #system size
    v0 = 0.3               #active velocity

    D = 0.00               #diffusion constant; T/drag
    Dr = 0.01              #rotational diffusion constant; T/drag_r

    epsilon = 1           #energy (an)isotropy; constant for now means isotropic
    beta = 2                #power in the potential; 2 for harmonic spring

    a_mor = 2.5
    D_mor = 0.08

    sigma_0 = 1
    aspect_ratio = 1.5      #ratio of major axis/minor axis; deformation only works if it is 1
    aspect_ratio_init = 1.5
    sigma_edge = sigma_0*aspect_ratio
    sigma_side = sigma_0
    lmda_major_0 = sigma_edge/2      #rest lambda of major axis
    lmda_minor_0 = sigma_side/2      #rest lambda of minor axis
    lmda_major_init = sigma_0*aspect_ratio_init/2
    lmda_minor_init = sigma_0/2
    R_0 = sigma_0

    tau = 1.0
    mu = 1.0
    K = 5.0

    # flag to switch on the deformations
    deformable = true

    area = pi*lmda_major_0*lmda_minor_0
    N = round(Int64,phi*L^2/area)
    #L =sqrt( N*area/phi)
    
    
    x = rand(Uniform(-L/2, L/2),N)
    y = rand(Uniform(-L/2, L/2),N)

    theta = rand(Uniform(0,2pi),N)
    println("number of particles : ",  N)


    lmda_major = ones(N)*lmda_major_init
    lmda_minor = ones(N)*lmda_minor_init
    sins = sin.(theta)
    coss = cos.(theta)
    #Lmda = np.array([[cos^2*lmda_major+sin^2*lmda_minor, sin*cos*(lmda_major-lmda_minor)], [sin*cos*(lmda_major-lmda_minor), sin^2*lmda_major+cos^2*lmda_minor]])
    i=1
    #display(  coss[i]^2*lmda_major[i]+sins[i]^2*lmda_minor[i])
    Lmbda = [  @MMatrix [coss[i]^2*lmda_major[i]+sins[i]^2*lmda_minor[i] sins[i]*coss[i]*(lmda_major[i]-lmda_minor[i]) ; sins[i]*coss[i]*(lmda_major[i]-lmda_minor[i]) sins[i]^2*lmda_major[i]+coss[i]^2*lmda_minor[i] ] for i=1:N ]


    #display(Lmbda)
    #Lmbda = hcat(@. [ coss[i]^2*lmda_major[i]+sins^2*lmda_minor[i] sins[i]*coss[i]*(lmda_major[i]-lmda_minor[i]); sins[i]*coss[i]*(lmda_major[i]-lmda_minor[i]) sins[i]^2*lmda_major[i]+coss[i]^2*lmda_minor[i]] for i=2:N])


    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3

    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    


    #dofevolvers = [inertial_evolver!]
    pair_forces =(Ellipse_2d_morse(1, D_mor, a_mor,1.),)
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_deformable_ellipse_evolver(1))
    global_dofevolvers = 
    field_dofevolvers = ()

    ϕ = phi
    ##
    # id
    # type

    # Shape of ellipse
    # first element represents major axis part

    # Lambda
    # Lambda0
    # sigma_0

    # K
    # mu
    # tau

    # v0
    # Dr
    
    # D

    # deformable

    # x
    # xuw
    # v
    # f
    # stress

    # p
    # q

    # zeta #Friction coefficient, only used icw overdamped integrators, first element consists for major axis
    # fact
    # fpas
    # ci
    ##
    initial_state = GBEllipses[ GBEllipses([i], [1], Lmbda[i], copy( Lmbda[i]), [sigma_0], [K], [mu], [tau], [v0], [Dr],[D], [deformable], [x[i], y[i],0],[x[i], y[i],0],[0, 0,0],[0, 0,0], [0 0; 0 0], [coss[i],sins[i],0], [0, 0,0], [1],[0, 0,0],[0, 0,0],[0, 0,0]) for i=1:N ];


    display(L)
    sizes = [L,L,2];
    print(sizes)
    initial_field_state= ()
    field_forces = ()
    field_updaters = ()

    #β=-1 interesting!
    external_forces = (ABP_perpendicular_angular_noise(1,[0,0,1]),ABP_3d_propulsion_force(1))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,4.);

    #Run integrationov
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plot}ting
    sim = Euler_integrator(system,dt, 10 ,fps=25,Tplot=1000,plot_functions=(plot_ellipses!,),plotdim=2, res=(1000,1000),record_folder_path="/Users/kammeraat/test_ellipses_tstop10/")#, plot_nematic_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim;

end

simulation() 

@profview sim = simulation()  