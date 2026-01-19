
include("AnalysisPipeline.jl")

using CairoMakie
CairoMakie.activate!()
raw_data_file = load_file("/Users/kammeraat/dwsa/single/simdata/v0_0.4/Dr_0.0/J_0.0/raw_data.h5")


figure_save_path= "/Users/kammeraat/dwsa/single/potential_plot/pot_plot.pdf"
begin
save_tax = raw_data_file["integration_info"]["save_tax"]

frame_numbers = 1:1#:length(save_tax)

frames = raw_data_file["frames"]
J = raw_data_file["system"]["forces"]["external"]["self_align_with_v_unit_force"]["β"]
Dr = frames["1"]["Dr"][1]
v0 = frames["1"]["v0"][1]   
t = Observable(0.)

field_C = frames["1"]["field_C"]
field_x_centers = frames["1"]["field_bin_centers_x"]
field_y_centers = frames["1"]["field_bin_centers_y"]


#Setup figure
f = Figure()
ax = Axis(f[1,1], aspect=DataAspect(),title = ("Surface plot potential"), xlabel="x", ylabel="y");
xlims!(ax, low=-2.5, high=2.5)
ylims!(ax, low=-1.5, high=1.5)

#potential 
Cmax = maximum(field_C)
Cmin = minimum(field_C)
surface!(ax,field_x_centers,field_y_centers,field_C.*0,colormap=Reverse(:gist_rainbow), rasterize = true, alpha=1.0, colorrange=(Cmin, Cmax), color=field_C,shading=NoShading)
Colorbar(f[1,2], limits = (Cmin, Cmax), label="Pot. energy U",colormap=Reverse(:gist_rainbow))

#scatter!(ax,points, color="black", alpha=0.2)
#meshscatter!(ax,x,y)

#directors
#Label(f[2,1],"System parameters: "*string(["$(key)=$(val)" for (key,val) in tag]), tellwidth=false, halign=:left, word_wrap = true)

display(f)
save("test.pdf", f, px_per_unit=100.0)
end
#save()
close(raw_data_file)



