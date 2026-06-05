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

struct GeneralField2d<:Field
    id::Int64
    type::Int64
    bin_centers::Vector{Vector{Float64}}
    C::Matrix{Float64}
    Cv::Matrix{Float64}
    Cf::Matrix{Float64}
    lbin::Float64
end

