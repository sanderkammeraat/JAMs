


include("AnalysisPipeline.jl")



file = load_file("/Users/kammeraat/sa_raw_data.h5")


frames = file["frames"]


Nint = length(frames["1"]["id"])
t = file["integration_info"]["save_tax"]


Nt = length(t)
x = zeros(Nint, Nt)
y = zeros(Nint, Nt)

vx = zeros(Nint, Nt)
vy = zeros(Nint, Nt)

px = zeros(Nint, Nt)
py = zeros(Nint, Nt)



@views for i in 1:Nt
    x[:,i] .= extract_frame_data_for_type("x", 1, frames[string(i)])
    y[:,i] .= extract_frame_data_for_type("y", 1, frames[string(i)])
    vx[:,i] .= extract_frame_data_for_type("vx", 1, frames[string(i)])
    vy[:,i] .= extract_frame_data_for_type("vy", 1, frames[string(i)])
    px[:,i] .= extract_frame_data_for_type("px", 1, frames[string(i)])
    py[:,i] .= extract_frame_data_for_type("py", 1, frames[string(i)])
end




begin
    f = Figure()
    ax =Axis(f[1,1], aspect=1)

    i=10
    scatter!(ax, x[:,i], y[:,i])

    display(f)
end

Lx = file["system"]["sizes"][1]
Ly = file["system"]["sizes"][2]


kxs = collect(range(2*pi/Lx, pi, step=2*pi/Lx))
kys = collect(range(2*pi/Ly, pi, step=2*pi/Ly))

#time index
i = 1

function fourier_transform_2d(x, y, kx, ky, vx)

    X = zeros(Complex,length(kx), length(ky))
    Kx = zeros(length(kx), length(ky))
    Ky =  zeros(length(kx), length(ky))

    for (i, kxi) in pairs(kx)
        for (j, kyj) in pairs(ky)

            Kx[i,j] = kxi
            Ky[i,j] = kyj

            for n in eachindex(vx)
                X[i,j] += vx[n] * exp(-1*im * kxi * x[n] -1*im * kyj * y[n] )
            end

        end
    end

    X2 = abs.(X) .* abs.(X)
    return Dict("Kx"=>Kx, "Ky"=>Ky, "kx"=>kx, "ky"=>ky, "X2"=>X2) 
end

i=1000
FT = fourier_transform_2d(x[i,:], y[i,:], kxs, kys, vx[i,:])


begin
    f = Figure()
    ax = Axis(f[1,1],)

    #heatmap!(ax, FT["kx"], FT["ky"], FT["X2"])
    X2 = FT["X2"]
    Kr =  sqrt.( FT["Kx"].^2 +  FT["Ky"].^2)
    kx =  FT["kx"]
    ky =  FT["ky"]
    scatter!(ax, [Kr[i,j] for i in eachindex(kx) for j in eachindex(ky)], [X2[i,j] for i in eachindex(kx) for j in eachindex(ky)])
    display(f)

end

begin
    f = Figure()
    ax = Axis(f[1,1],aspect=1)

    heatmap!(ax, FT["kx"], FT["ky"], FT["X2"])

    display(f)

end




close(file)