abstract type Field end

struct FuelField2d<:Field
    x_bin_centers::Vector{Float64}
    y_bin_centers::Vector{Float64}
    C::Array{Float64}

    Cv::Array{Float64}
    Cf::Array{Float64}
end

