

#Run this script with one thread from terminal, n will run the different parameter space points over n cores
using ArgParse
using Distributed
simargs = ArgParseSettings()

@add_arg_table simargs begin
    "--ncores", "-n"
        arg_type=Int64

end

parsed_params = parse_args(simargs)
n=parsed_params["ncores"]

addprocs(n)


@everywhere include(joinpath("..","src","Engine.jl"))

@everywhere function simulation(J, Dr,save_folder_path)

    external_forces = ( ABP_3d_propulsion_force(1), self_align_with_v_unit_force(1,J),ABP_perpendicular_angular_noise(1,[0,0,1]))

    pair_forces = [soft_disk_force([1, 2],[1. 1.; 1. 1.])]

    local_dofevolvers = [overdamped_xvf_evolver(1), overdamped_pq_evolver(1)]
    global_dofevolvers = []
    field_dofevolvers = []


    #First make stair
    
    Nlin=4
    Nrows = 2*Nlin
    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]

    xs = []
    ys = []
    r=1.0
    l = 2# 2*r

    typess = []

    for row in 1:Nlin

        push!(xs, [xi-(Nlin+1)*l/2-(row-1)*l/2 for xi in l.*range(1,Nlin+row-1) ])

        push!(ys, [-(Nlin-row)/2*l*sqrt(3) for n in xs[row] ])

        if row==1
            rowtypes=2*ones(Nlin)

        else
            rowtypes = ones(length(xs[row]))
            rowtypes[1]=2
            rowtypes[end]=2

        end
        push!(typess,rowtypes )

        
    end

    for row in 1:Nlin-1

        push!(xs, xs[Nlin-row])
        push!(ys, ys[Nlin][1].-ys[Nlin-row])

        push!(typess, typess[Nlin-row])

    end


    x = vcat(xs...)
    y = vcat(ys...)
    types= vcat(typess...)
    N = length(x)


    poly=0.15e-2
    Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)

    initial_state = Union{PolarParticle3d,ConfinedPolarParticle3d}[]
    id=1
    #First loop over normal particles
    for i=1:N

        if types[i]==1
            
            push!(initial_state, PolarParticle3d([id], [types[i]], [1], [1], [Rs[i]], [0.01], [Dr], [x[i] , y[i],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]))
            id+=1
        end
    end
    for i=1:N

        if types[i]==2
            push!(initial_state,ConfinedPolarParticle3d([id],[2], [1],[1], [Rs[i]], [0], [0.01], [x[i] , y[i],0],[0.,0.,0.],[0,0,0], [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]))
            id+=1
        end
    end


    size = [Nrows*l+2*l,Nrows*l+2*l,1];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(size, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, local_dofevolvers,global_dofevolvers, field_dofevolvers, false,2.5);

    #Run integration
    sim = Euler_integrator(system,1e-2, 5e3,Tsave=100,Tplot=100, save_functions = [save_2d_polar_p!],save_folder_path=save_folder_path, fps=120, plot_functions=(plot_disks_orientation!,plot_directors!, plot_velocity_vectors!), plotdim=2); 
    return sim

end




Drs = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.]
J=1

@sync @distributed for Dr in Drs

    display("Running")

    save_folder_path = joinpath("/","data1","kammeraat","sa","varyDr","J_$J", "Dr_$Dr");
    print(save_folder_path)

    sim = simulation(J,Dr, save_folder_path);
end
