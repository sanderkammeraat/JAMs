abstract type Field end

#Include an id to the field to have its own rng.

struct FuelField2d<:Field
    id::Int64
    type::Int64
    bin_centers::Vector{Vector{Float64}}
    C::Matrix{Float64}
    Cv::Matrix{Float64}
    Cf::Matrix{Float64}
end

struct GeneralField2d{Tbc}<:Field
    id::Int64
    type::Int64
    bin_centers::Tbc
    C::Matrix{Float64}
    Cv::Matrix{Float64}
    Cf::Matrix{Float64}
    lbin::Float64
end

#Have cpu viewable matrix and a gpu viewable matrix
struct GPUField2d{Tc, Tcv, Tcf}<:Field
    id::Int64
    type::Int64
    bin_centers::Vector{Vector{Float64}}
    C::Matrix{Float32}
    Cv::Matrix{Float32}
    Cf::Matrix{Float32}

    C_GPU::Tc
    Cv_GPU::Tcv
    Cf_GPU::Tcf

    lbin::Float32
end

