
function extract_data_for_type(datakey, type, frame)

    return frame[datakey][ frame["type"].==type ]

end


function construct_n_ij_projector(i,j,x,y)


    r_ij_0_vec = [x[j]-x[i], y[j] - y[i]]

    r_ij_0_norm = norm(r_ij_0_vec)

    n_ij =r_ij_0_vec./r_ij_0_norm

    n_ij_projector = [ n_ij[1]*n_ij[1] n_ij[1]*n_ij[2] ; n_ij[2]*n_ij[1] n_ij[2]*n_ij[2]]

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

function construct_M_ij(i,j,x,y, k, R,type)

    M_ij= zeros(2,2)

    if i!=j
        n_ij_projector, r_ij_0_norm = construct_n_ij_projector(i,j,x,y)

        M_ij.+= u_ij_1(i,j,r_ij_0_norm, k, R, type) / r_ij_0_norm .* (I - n_ij_projector)

        M_ij.+= u_ij_2(i,j,r_ij_0_norm, k, R, type) .* n_ij_projector

    end

    return M_ij

end



function construct_D(x0,y0, k, R, type)

    M=zeros(2*length(x0), 2*length(x0))

    for i in eachindex(x0)

        for j in i:length(x0)

            M[2i-1:2i,2j-1:2j] = construct_M_ij(i, j, x0, y0, k , R,type)
        end

    end
    D = -(M + transpose(M))

    for i in eachindex(x0)

        for j in eachindex(x0)

            #In principle unnecessary because Mii =0, but just to be sure
            if i!=j
                D[2i-1:2i,2i-1:2i]-= D[2i-1:2i,2j-1:2j]
            end
        end
    end
    return Symmetric(D)
end