include(joinpath("..","src","Engine.jl"))
using Random, Distributions

#unphysical system 
function simulation()
    lx = 4

    ly=lx

   

    lxx = lx*1
    lyy=lxx

    mx = 10*lx
    my=mx

    x = -mx/2:lxx:mx/2
    y = -my/2:lyy:my/2

    Nx = length(x)
    Ny = length(y)
    N = Nx*Ny
    initial_state= PolarParticle3d[]
    id=0
    for i in 1:Nx
        for j in 1:Ny
            
            id+=1

            push!(initial_state,PolarParticle3d([id],[1],[1],[1],[1],[0.01],[0.01],[x[i],y[j],0],[0.,0.,0.],[0,0,0],[0,0,0],[0,0,0],normalize([rand(Normal(0, 1)),rand(Normal(0, 1)),0]),[0,0,0],[0,0,0]) )

        end


    end

    external_forces = (ABP_3d_propulsion_force(1), self_align_with_v_force(1,1),ABP_perpendicular_angular_noise(1,[0,0,1]))

    dofevolvers=[overdamped_evolver!]

    k_network = zeros(N, N)
    k=1e-1
    for (i, p_i) in pairs(initial_state)

        for (j,p_j) in pairs(initial_state)
            if (abs(p_j.x[1]-p_i.x[1])==lxx) && (abs(p_i.x[2]-p_j.x[2])==0)

                k_network[i,j] = k

            end

            if (abs(p_i.x[2]-p_j.x[2])==lyy) && (abs(p_i.x[1]-p_j.x[1])==0)

                k_network[i,j] = k

            end


        end

    end


    #display(k_network)

    pair_forces = [spring_network_2d_force(1,lx, k_network)]

    sizes = [3*findmax(x)[1],3*findmax(y)[1],1.];
    initial_field_state=[]
    field_forces = []
    field_updaters = []


    system = System(sizes, initial_state,initial_field_state, external_forces, pair_forces,field_forces, field_updaters, dofevolvers, true,2.);

    sim = Euler_integrator(system,1e-2, 1e5, Tplot=5e1, fps=120, plot_functions=(plot_disks_orientation!, plot_directors!, plot_velocity_vectors!), plotdim=2); 
end


sim = simulation()