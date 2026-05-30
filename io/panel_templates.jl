inch = 96
pt = 4/3
cm = inch / 2.54
fontsize = 7pt

draft = Theme(
    fontsize = fontsize,
    figure_padding = 10,
    font = "Helvetica",
    fonts = (; 
        regular = "Helvetica", 
        bold = "Helvetica Bold", 
        italic = "Helvetica Oblique", 
        bold_italic = "Helvetica Bold Oblique",
        math = "Helvetica"
    ),
    Axis = (
        xlabelsize = fontsize,  
        xticklabelsize = fontsize,
        yticklabelsize = fontsize,  
        ylabelsize = fontsize,
        xtickalign = 1,        
        ytickalign = 1,
        xticksmirrored = true,
        yticksmirrored = true,
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminortickalign = 1,
        yminortickalign = 1,
        xminorticks = IntervalsBetween(5),
        yminorticks = IntervalsBetween(5),
        width= 4cm,
        height= 3cm,
        xgridvisible = false, 
        ygridvisible = false,
    ),
    Legend = (
        labelsize = fontsize,
        patchsize = (10, 10),
        framevisible = false,
        padding = (2, 2, 2, 2),
        tellheight=false,
        tellwidth=true,
        halign=:left
    )
)