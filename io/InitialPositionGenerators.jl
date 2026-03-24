
using Distributions

function Random_in_hexagon(N, xb, yb, Rb)

    #Get maximum radius of boundary particles in case of polydispersity
    Rbmax = maximum(Rb)

    xbmin = minimum(xb)+ 2* Rbmax
    xbmax = maximum(xb) - 2*Rbmax

    ybmin = minimum(yb)+ 2*Rbmax
    ybmax = maximum(yb) - 2*Rbmax


    x = zeros(N)
    y = zeros(N)

    for i in 1:N

        #Trial point can be limited to bounding box, to reduce number of trials needed


        #Then reject any points inside the outer triangles complementing the hexagon in the bounding box

        pass = false

        while pass == false

            x_trial = rand(Uniform(xbmin, xbmax ))

            y_trial = rand(Uniform(ybmin, ybmax ))

            #Let's give it a shot
            pass = true

            #top left corner:
            pass *= ( y_trial  < sqrt(3) * (x_trial -xbmin))


            #top right corner
            pass *= y_trial  < - sqrt(3) * (x_trial - xbmax)
            
            # bottom left corner 

            pass *= y_trial  - ybmin > - sqrt(3) * (x_trial - xbmin) + ybmax

            # bottom right corner 
            pass *= y_trial  - ybmin >  sqrt(3) * (x_trial - xbmax) + ybmax

            if pass
                x[i] = x_trial
                y[i] = y_trial
            end

        end


    end

    return x, y


end


#Testing
# begin
# Nlin=50
# Nrows = 2*Nlin

# xs = []
# ys = []
# r=1.0
# ϕ=1.1
# l = 2* sqrt(pi*sqrt(3)/(6*ϕ))# 2*r


# typess = []

# for row in 1:Nlin

#     push!(xs, [xi-(Nlin+1)*l/2-(row-1)*l/2 for xi in l.*range(1,Nlin+row-1) ])

#     push!(ys, [-(Nlin-row)/2*l*sqrt(3) for n in xs[row] ])

#     if row==1
#         rowtypes=2*ones(Nlin)

#     else
#         rowtypes = ones(length(xs[row]))
#         rowtypes[1]=2
#         rowtypes[end]=2

#     end
#     push!(typess,rowtypes )

    
# end

# for row in 1:Nlin-1
#     push!(xs, xs[Nlin-row])
#     push!(ys, ys[Nlin][1].-ys[Nlin-row])

#     push!(typess, typess[Nlin-row])

# end

# const x = vcat(xs...)
# const y = vcat(ys...)
# types= vcat(typess...)
# const N = length(x)
# println(N)

# poly=0.15*1e0
# const Rs = rand(Uniform((1-poly)*r, (1+poly)*r),N)
# end


# xi, yi = Random_in_hexagon(100000, x, y, Rs)

# using GLMakie

# GLMakie.activate!()

# scatter(x, y)
# scatter!(xi, yi)

#Thanks to Gabriel Martin
function stacked_polymers_at_angle(N_in_pol, R, pf, f_eq_stretch_force, L_0; tilt_angle = nothing, random_polarity = false)

    # initialization

    S = pi * R^2

    S_overlap = 2*(R*R*acos(f_eq_stretch_force*R) - (f_eq_stretch_force*R)^2*tan(acos(f_eq_stretch_force*R)))

    Spol = S*N_in_pol-S_overlap*(N_in_pol-1)

    Npols = convert(Int64, floor(L_0^2 * pf / Spol))

    L = sqrt(Npols*Spol/pf)

    M = Npols * N_in_pol

    x = zeros(M)
    y = zeros(M)
    radii = R * ones(M)

    pol_ids = zeros(M)

    ids_in_pol= zeros(M)

    Lpols = 2*R + (N_in_pol - 1) * 2*R * f_eq_stretch_force

    if isnothing(tilt_angle)
        tilt_angle = atan(L/Npols/Lpols)   # The angle along which to place the particles to maximize the spread of particles on the torque
    end
    
    id = 1

    for (i, pol_index) in enumerate(shuffle!(collect(1:Npols)))

        flip = random_polarity ? rand((true,false)) : false   #based on the random_polarity bool, if true choose random, if false use false

        for j in 1:N_in_pol
            
            pol_ids[id] = pol_index

            ids_in_pol[id] = j

            index = j - 1

            if flip
                index = N_in_pol - j
            end

            x[id] = ((i-1)*Lpols + index*2*R*f_eq_stretch_force)*cos(tilt_angle) % L - L/2
            y[id] = ((i-1)*Lpols + index*2*R*f_eq_stretch_force)*sin(tilt_angle) % L - L/2

            id+=1

        end

    end

    x[1,1] *= 0.99
    y[1,1] *= 0.99

    return x, y, radii, pol_ids, ids_in_pol, L, Npols

end