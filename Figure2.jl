
alpha_trans=1
anno2 = Dict("Podo" => ("magenta1",alpha_trans), 
                "HealthyPT"=>("green3",alpha_trans),
                "InjPT"=>("#f92874",alpha_trans),
                "TAL"=>("lightslateblue",alpha_trans),
                "DCT"=>("blue",alpha_trans),
                "CD-PC"=>("turquoise1",alpha_trans),
                "CD-IC"=>("#924cfa",alpha_trans),
                "vEC"=>("firebrick",alpha_trans),
                "gEC"=>("dodgerblue",alpha_trans),
                "Fib"=>("#edff4d",alpha_trans),
                "JGA"=>("sienna2",alpha_trans),
                "Immune"=>("darkgreen",alpha_trans),
                "Uro"=>("black",alpha_trans));
p3=spd.plot_cell_polygons(sham, "celltype"; 
    anno_color=anno2,x_lims=(0,35000), 
    y_lims=(0,40000),canvas_size=(5000,6000),
    stroke_color="gray80")
