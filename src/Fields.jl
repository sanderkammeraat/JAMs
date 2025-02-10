abstract type Field end

struct FuelField2d<:Field
    type::Int64
    bin_centers::Vector{Vector{Float64}}
    C::Array{Float64}
    Cv::Array{Float64}
    Cf::Array{Float64}
end

struct GeneralField2d<:Field
    type::Int64
    bin_centers::Vector{Vector{Float64}}
    C::Array{Float64}
    Cv::Array{Float64}
    Cf::Array{Float64}
end

