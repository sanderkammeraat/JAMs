include(joinpath("..","src","Engine.jl"))

function simulation(v0, Dr, J, Tplot, Tsave)

    

    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,1.2,0.3,0.3,true), pair_nematic_alignment_force(1,2.5,0.1))
    #type, torque, rfact, kpar, kper
    #1.3 -1 0 0.3
    #pair_forces = (soft_disk_force(1,1),pairAN_force(1,true,1.3, 1, 0., 0.3), pair_nematic_alignment_force(1,2.5,0.15))
    pair_forces =[]#(soft_disk_force(1,0),)


    #dofevolvers = [inertial_evolver!]
    local_dofevolvers = (overdamped_xvf_evolver(1),overdamped_pq_xyc_evolver(1))
    global_dofevolvers = []
    field_dofevolvers = []

    N=1000
    ϕ = 1.0
    poly=15e-6
    Rs = rand(Uniform(1-poly, 1+poly),N)

    L =  5
    Lx = L
    Ly = L
    #v0 = 0.2
    #Dr = 0.0
    sigma = 1
    initial_state = PolarParticle3d[ PolarParticle3d([i],[1], [1], [1], [Rs[i]], [v0], [Dr], [-1+sigma*rand(Normal(0, 1)) , sigma*rand(Normal(0, 1)),0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) for i=1:N ];
    lbin=0.1
    z_bin_centers = [0.]
    x_bin_centers = []
    x_bin_centers = append!(x_bin_centers,range(start=-2*Lx/3, stop=2*Lx/3+0.1*lbin, step=lbin).+lbin/2)

    y_bin_centers = []
    y_bin_centers = append!(y_bin_centers,range(start=-2*Ly/3, stop=2*Ly/3+0.1*lbin, step=lbin).+lbin/2)
    bin_centers = [x_bin_centers, y_bin_centers,z_bin_centers]

    C = zeros(length(x_bin_centers), length(y_bin_centers))

    xa = -1
    xb = 1
    ka = 2
    kb = 2
    for (i, xi) in pairs(x_bin_centers)

        for (j, yj) in pairs(y_bin_centers)

            C[i,j] = -exp( - ka/2* ( (xi - xa)^2 + yj^2)) - exp( - kb/2* ((xi - xb)^2 +yj^2))
        end
    end


    initial_field_state=[FuelField2d(1,1,bin_centers,C, C.*0, C.*0)]

    display(L)
    sizes = [L,L,4];
    #initial_field_state=[]
    field_forces = []
    field_updaters = []

    #β=-1 interesting!
    external_forces = (external_double_gaussian_force(1,ka,kb, [xa,0,0], [xb,0,0]),ABP_perpendicular_angular_noise(1,[0,0,1]),self_align_with_v_unit_force(1,J),ABP_3d_propulsion_force(1))

    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, true,3.);

    #Run integrationov
    #Use plot_disks! for nice visualss
    #Use plot_points! for fast plotting

    #save_folder = "/Users/kammeraat/dwsa/single/simdata/v0_$v0/Dr_$Dr/J_$J/"
    save_folder = "/Users/kammeraat/mounting/data2_kammeraat/test_saving/single/simdata/v0_$v0/Dr_$Dr/J_$J/"
    sim = Euler_integrator(system,0.01,1, Tsave=Tsave, fps=60,Tplot=Tplot,plot_functions=(plot_potential!,plot_points!, plot_velocity_vectors!), plotdim=2, save_folder_path = save_folder, save_functions = (save_2d_polar_p!,save_single_2d_field!)); 
    return sim;

end

#v0, Dr, J
simulation(0.6, 0.01 ,1, 1, 1)

Js = [0, 0.1, 0.5, 1.0, 2.0]

Drs = [0.0, 0.01, 0.1]

for Dr in Drs
    for J in Js
        simulation(0.4, Dr, J, nothing, 1)
        #sleep(2)
    end
end
