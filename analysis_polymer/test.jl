using CairoMakie

f = Figure()
ax = Axis(f[1,1], title="ploup")

fon(x) = 2*x


t = [1,2,3,4,5,6,7,8,9]
lines!(t, fon(t))
save("C:\\Users\\gabri\\Documents\\plapla.pdf", f)