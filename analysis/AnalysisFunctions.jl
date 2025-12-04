
using FFTW


function temporal_Fourier_transform(dt, x;  min_t_ind=1, output_not_avg=false)

    #Calculate the Fourier transform along the time axis for a real signal

    @views xf = x[:, min_t_ind:end ]
    Xf =rfft(xf, 2 )

    #Times 2 π to get angular frequency
    w = convert(Array, rfftfreq( size(xf)[2], 1/dt).*2*pi)


    Xf2 = abs.(Xf).^2

    @views pavg_Xf2 = mean(Xf2, dims=1)[1,:]

    @views pstd_Xf2 = std(Xf2, dims=1)[1,:]

    #Standard error to the mean
    @views pste_Xf2 = pstd_Xf2 ./ sqrt(size(Xf2)[1])




    if output_not_avg

            maxval, maxind = findmax(Xf2, dims=2)
            w_max = [w[ci[2]] for ci in maxind[:,1] ]

        FT = Dict("Xf2"=>Xf2,  "w"=>w,"max_X2"=>maxval, "max_X2_ind" =>maxind, "w_max"=>w_max )

    else
        maxval, maxind = findmax(pavg_Xf2)
        w_max = w[maxind]
        FT = Dict("pavg_X2"=>pavg_Xf2, "pstd_X2"=>pstd_Xf2, "pste_X2"=>pste_Xf2, "w"=>w, "max_X2"=>maxval, "max_X2_ind" =>maxind, "w_max"=>w_max )
    end
    return FT

end

function secondary_temporal_Fourier_transform(dt, C)


        #Calculate the Fourier transform along the time axis for a real signal
        Xf =rfft(C)
    
        #Times 2 π to get angular frequency
        w = rfftfreq( length(C), 1/dt).*2*pi
    
    
        Xf2 = abs.(Xf).^2
    
    
        maxval, maxind = findmax(Xf2)
        w_max = w[maxind]
    
    
        FT = Dict("X2"=>Xf2, "w"=>w, "max_X2"=>maxval, "max_X2_ind" =>maxind, "w_max"=>w_max )
    
        return FT

end


function spatial_p_correlation(binsize, maxbin_center, x,y, px, py; min_t_ind=1, max_t_ind=nothing, every_n=nothing)


    rbin_edges = prepend!(collect(range(start=0, step=binsize, stop=maxbin_center)),[0])

    rbin_edges2 = rbin_edges.^2

    rbin_centers = (rbin_edges[2:end] + rbin_edges[1:end-1])/2

    Nbin = length(rbin_centers)

    Nt = size(px)[2]

    max_t_ind_set = isnothing(max_t_ind) ? Nt : max_t_ind

    every_n_set = isnothing(every_n) ? 1 : every_n

    C = ones(Nt, Nbin)*NaN

    counts = zeros(Nt, Nbin)

    @showprogress dt = 1 desc="spatial correlation" showspeed=true for i in min_t_ind:every_n_set:max_t_ind_set
        for p1 in 1:size(px)[1]

            for p2 in p1:size(px)[1]

                Δr2 = (x[p1, i]- x[p2,i])^2 +   (y[p1, i]- y[p2,i])^2

                for bin in eachindex(rbin_centers)

                    if (Δr2<= rbin_edges2[bin+1] && Δr2> rbin_edges2[bin]) || (p1==p2 && bin==1)

                        if isnan(C[i, bin])
                            C[i, bin]=0
                            counts[i, bin]=0
                        end

                        C[i, bin] += px[p1, i]* px[p2, i] + py[p1, i] * py[p2, i]
                        counts[i,bin] += 1
                    end

                end

            end
        end
    end

    C.=C./counts
    return Dict("rbc"=>rbin_centers,"rbe"=>rbin_edges,"C"=> C)

end



function spatiotemporal_p_correlation(binsize, maxbin_center, x0,y0, px, py; min_t_ind=1, max_t_ind=nothing, every_n=nothing)

    rbin_edges = prepend!(collect(range(start=0, step=binsize, stop=maxbin_center)),[0])

    rbin_edges2 = rbin_edges.^2

    rbin_centers = (rbin_edges[2:end] + rbin_edges[1:end-1])/2

    Nbin = length(rbin_centers)

    Nt = size(px)[2]

    C = ones(Nt, Nbin)*NaN

    counts = zeros(Nt, Nbin)

    max_t_ind_set = isnothing(max_t_ind) ? Nt : max_t_ind

    #every_n_set = isnothing(every_n) ? 1 : every_n

    #particle 1 loop
    @showprogress dt = 1 desc="spatiotemporal correlation" showspeed=true for p1 in 1:size(px)[1]

        #particle 2 loop
        for p2 in 1:size(px)[1]

            Δr2 = (x0[p1]- x0[p2])^2 +   (y0[p1]- y0[p2])^2

            #Loop over spatial bin
            for bin in eachindex(rbin_centers)

                if (Δr2<= rbin_edges2[bin+1] && Δr2> rbin_edges2[bin]) || (p1==p2 && bin==1)

                    for i in min_t_ind:max_t_ind_set

                        dij = abs(i)

                        if isnan(C[dij, bin])
                            C[dij, bin]=0
                            counts[dij, bin]=0
                        end

                        C[dij, bin] += px[p1, min_t_ind]* px[p2, i] + py[p1, min_t_ind] * py[p2, i]
                        counts[dij,bin] += 1
                    end
                end

            end

        end
    end
    C.=C./counts
    return Dict("rbc"=>rbin_centers,"rbe"=>rbin_edges,"C"=> C)
end

# function auto_correlation(t, vx, vy)

#     Nt = length(t)
#     Δt = zeros(Nt)
#     C = zeros(Nt, Nt)
#     N_in_bin = zeros(Nt, Nt)

#     @showprogress desc="Autocorrelation" dt=1 showspeed=true for i in 1:Nt-1

#         for j in i:Nt
#             C[i,j-i+1] = mean(px[:,j].* px[:,i] .+ py[:,j].* py[:,i])






function auto_correlation(t, px, py; normalized=false, minrow=1, maxrow=nothing)


    Nt = length(t)

    Δt = zeros(Nt)

    C = zeros(Nt, Nt).*NaN

    @showprogress desc="Autocorrelation" dt=1 showspeed=true for i in 1:Nt-1

         @Threads.threads for j in i:Nt
            @views C[i,j-i+1] = mean(px[:,j].* px[:,i] .+ py[:,j].* py[:,i])
            Δt[j-i+1] =  t[j] - t[i]
        end
        

    end

    Cavg = zeros(Nt-minrow+1)
    for j in 1:Nt-minrow+1
        if isnothing(maxrow)
            Cavg[j] = mean( filter!( e->!isnan(e) ,C[minrow:end,j] ) )
        else
            Cavg[j] = mean( filter!( e->!isnan(e) ,C[minrow:maxrow,j] ) )
        end
    end

    return Dict("Cavg"=>Cavg, "t"=>t, "deltat"=> Δt)
end


function auto_correlation_v2(t, px, py; minrow=1)

    Cp = zeros(size(px[:,minrow:end]))
    delta_t = zeros(length(t[minrow:end]))

    for i in 1:size(Cp)[1]
        for j in 1:size(Cp)[2]

            Cp[i,j] = px[i,minrow - 1 +j ] * px[i, minrow] + py[i,minrow - 1 +j ] * py[i, minrow]
            delta_t[j] = t[minrow-1+j] - t[minrow]
        end
    end
    Cavg = mean(Cp, dims=1)[1,:]

    return Dict("Cavg"=>Cavg, "t"=>t, "deltat"=>delta_t)


end
## Dynamical matrix analysis
# Helper functions
function construct_n_ij_projector!(n_ij_projector, i,j,x,y)


    r_ij_0_vec = @MVector [x[j]-x[i], y[j] - y[i]]

    r_ij_0_norm = norm(r_ij_0_vec)

    n_ij =r_ij_0_vec./r_ij_0_norm

    n_ij_projector[1,1] = n_ij[1]*n_ij[1]
    n_ij_projector[1,2] = n_ij[1]*n_ij[2]
    n_ij_projector[2,1] = n_ij[1]*n_ij[2]
    n_ij_projector[2,2] = n_ij[2]*n_ij[2]

    return n_ij_projector, r_ij_0_norm
end

function u_ij_1(i,j,r_ij_0_norm, k, R,type)
    #If not the same
    if i==j
        return 0.

    #If different
    elseif i!=j
        #If in 
        a = R[i]+R[j]
        if r_ij_0_norm<=a
            return k[type[i],type[j]] * (r_ij_0_norm - a)

        else
            return 0.
        end
    end
end

function u_ij_2(i,j,r_ij_0_norm, k, R,type)
    #If the same
    if i==j
        return 0.

    #If different
    elseif i!=j
        #If in 
        a = R[i]+R[j]
        if r_ij_0_norm<=a
            return k[type[i],type[j]]
        else
            return 0.
        end
    end
end

@views function construct_M_ij!(M_ij,n_ij_projector, i,j,x,y, k, R,type)

    M_ij .*=  0

    if i!=j
        n_ij_projector, r_ij_0_norm = construct_n_ij_projector!(n_ij_projector,i,j,x,y)

        M_ij.+= u_ij_1(i,j,r_ij_0_norm, k, R, type) / r_ij_0_norm .* (I - n_ij_projector)

        M_ij.+= u_ij_2(i,j,r_ij_0_norm, k, R, type) .* n_ij_projector

    end

    return M_ij

end


# Actual D construction
@views function construct_D(x0,y0, k, R, type)

    M=zeros(2*length(x0), 2*length(x0))

    D = copy(M)

    #Allocate once
    M_ij_0=  @MMatrix zeros(2,2)

    n_ij_projector =  @MMatrix zeros(2,2)
    

    @showprogress for i in eachindex(x0)

        for j in i:length(x0)

            M[2i-1:2i,2j-1:2j] .= construct_M_ij!(M_ij_0, n_ij_projector,i, j, x0, y0, k , R,type)
        end

    end
    D  .= -(M + transpose(M))

    @showprogress for i in eachindex(x0)

        for j in eachindex(x0)

            #In principle unnecessary because Mii =0, but just to be sure
            if i!=j
                 D[2i-1:2i,2i-1:2i] .-= D[2i-1:2i,2j-1:2j]
            end
        end
    end
    return Symmetric(D)
end






function diagonalize_D(D)

    eigenfact = eigen(D)
    return Dict("eigvals"=>eigenfact.values, "eigvecs"=>eigenfact.vectors)

end





@views function project_on_eigvecs(eigvecs, xinterior, yinterior)

    Neigvecs = size(eigvecs)[2]
    Nt = size(xinterior)[2]

    projs = zeros(Neigvecs, Nt)

    @showprogress dt = 1 desc="Projection on eigvecs" showspeed=true   for i in 1:Nt

         Threads.@threads for j in 1:Neigvecs
             projs[j,i] =  project_on_eigvec(eigvecs[:,j], xinterior[:,i], yinterior[:,i])

        end


    end

    return projs
end

@views function project_on_eigvec(eigvec, xi, yi)

    xy_zip = collect(Iterators.flatten(zip(xi, yi)))
    proj = sum( xy_zip .* eigvec)
    return proj
end


function unwrap(angles)
    θpc = zeros(size(angles))
    θpc[:,1] = angles[:,1]

    for p in 1:size(angles)[1]

        θp = angles[p,:]

        d = diff(θp)

        for (i, di) in pairs(d) 
            if -pi<di<pi
                θpc[p,i+1] =  θpc[p,i] + di

            elseif di>pi
                θpc[p,i+1] =  θpc[p,i] + di-2*pi

            elseif di<-pi
                θpc[p,i+1] =  θpc[p,i] + di+2*pi
            end
        end
    end
    return θpc
end