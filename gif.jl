using Plots
using Plots.PlotMeasures
using DelimitedFiles

cs = readdlm("cmake-build-debug/90%HR./test/dane_171.5_0.1185_80_350_1.37_1_0.073.txt")

plot(cs[:,1], cs[:,2],
    ylabel = "Cross-section [pb]",
    xlabel = "E [GeV]",
    size = (1200, 1000),
    tickfontsize = 16,
    guidefontsize = 20,
    legendfont = 16,
    linewidth = 5,
    bottom_margin = 5mm,
    legend = nothing,
    left_margin = 5mm,
    framestyle = :box,
    color = :red)

d = Dict(cs[:,1] .=> cs[:,2])

f = open("cmake-build-debug/ewolucja.txt", "r")
anim = @animate for i in 1:35
     # read a new / next line for every iteration
     plot(cs[:,1], cs[:,2],
         ylabel = "Cross-section [pb]",
         xlabel = "E [GeV]",
         size = (1200, 1000),
         tickfontsize = 16,
         guidefontsize = 20,
         legendfont = 16,
         linewidth = 5,
         bottom_margin = 5mm,
         legend = nothing,
         left_margin = 5mm,
         framestyle = :box,
         color = :red)

     s = parse.(Float64, split(readline(f)))
     c = zeros(Float64, length(s))
     for i in 1:length(s)
         c[i] = get(d, s[i], NaN)
     end
     scatter!(s, c,
              markersize = 10)
end

gif(anim, "ewolucja.gif", fps=4)
