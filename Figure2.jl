import CellScopes as cs
## whole kidney
alpha_trans=1
anno2 = Dict("Podo" => ("magenta1",alpha_trans), "HealthyPT"=>("green3",alpha_trans), "InjPT"=>("#f92874",alpha_trans),"TAL"=>("lightslateblue",alpha_trans),
            "DCT"=>("blue",alpha_trans),"CD-PC"=>("turquoise1",alpha_trans),"CD-IC"=>("#924cfa",alpha_trans),"vEC"=>("firebrick",alpha_trans),
             "gEC"=>("dodgerblue",alpha_trans),"Fib"=>("#edff4d",alpha_trans),"JGA"=>("sienna2",alpha_trans),"Immune"=>("darkgreen",alpha_trans), "Uro"=>("black",alpha_trans));
cs.plot_cell_polygons(sham, "celltype"; 
    cell_colors=anno2,x_lims=(0,35000), 
    y_lims=(0,40000), width=5000,height=6000,
    stroke_color="gray80")

## artery
anno2 = Dict("Podo" => ("white",alpha_trans), 
                "HealthyPT"=>("white",alpha_trans),
                "InjPT"=>("white",alpha_trans),
                "TAL"=>("white",alpha_trans),
                "DCT"=>("white",alpha_trans),
                "CD-PC"=>("white",alpha_trans),
                "CD-IC"=>("white",alpha_trans),
                "vEC"=>("firebrick",alpha_trans),
                "gEC"=>("white",alpha_trans),
                "Fib"=>("#edff4d",alpha_trans),
                "JGA"=>("white",alpha_trans),
                "Immune"=>("white",alpha_trans),
                "Uro"=>("white",alpha_trans));
cs.plot_cell_polygons(sham, "celltype"; 
    cell_colors=anno2, x_lims=(21900,24500), y_lims=(11500,19500),
    width=300,height=900,stroke_color="gray80")

## podocyte
anno2 = Dict("Podo" => ("magenta1",alpha_trans), 
                "HealthyPT"=>("white",alpha_trans),
                "InjPT"=>("white",alpha_trans),
                "TAL"=>("white",alpha_trans),
                "DCT"=>("white",alpha_trans),
                "CD-PC"=>("white",alpha_trans),
                "CD-IC"=>("white",alpha_trans),
                "vEC"=>("white",alpha_trans),
                "gEC"=>("dodgerblue",alpha_trans),
                "Fib"=>("white",alpha_trans),
                "JGA"=>("sienna2",alpha_trans),
                "Immune"=>("white",alpha_trans),
                "Uro"=>("white",alpha_trans))
cs.plot_cell_polygons(sham, "celltype"; 
    cell_colors=anno2, x_lims=(17800,20000), y_lims=(8700,9600), 
    width = 600,height = 250,stroke_color="gray80")

## cortex
anno2 = Dict("Podo" => ("magenta1",alpha_trans), 
                "HealthyPT"=>("green3",alpha_trans),
                "InjPT"=>("#f92874",alpha_trans),
                "TAL"=>("lightslateblue",alpha_trans),
                "DCT"=>("blue",alpha_trans),
                "CD-PC"=>("turquoise1",alpha_trans),
                "CD-IC"=>("#924cfa",alpha_trans),
                "vEC"=>("firebrick",alpha_trans),
                "gEC"=>("dodgerblue",alpha_trans),
                "Fib"=>("#edff4d",0.5),
                "JGA"=>("sienna2",alpha_trans),
                "Immune"=>("darkgreen",alpha_trans),
                "Uro"=>("black",alpha_trans))
cs.plot_cell_polygons(sham, "celltype"; 
    cell_colors=anno2, x_lims=(23000,24800), y_lims=(7400,9000), 
    width=450,height=400,stroke_color="gray80")

## outer medulla
anno2 = Dict("Podo" => ("magenta1",alpha_trans), 
                "HealthyPT"=>("green3",alpha_trans),
                "InjPT"=>("#f92874",alpha_trans),
                "TAL"=>("lightslateblue",alpha_trans),
                "DCT"=>("blue",alpha_trans),
                "CD-PC"=>("turquoise1",alpha_trans),
                "CD-IC"=>("#924cfa",alpha_trans),
                "vEC"=>("firebrick",alpha_trans),
                "gEC"=>("dodgerblue",alpha_trans),
                "Fib"=>("#edff4d",0.5),
                "JGA"=>("sienna2",alpha_trans),
                "Immune"=>("darkgreen",alpha_trans),
                "Uro"=>("black",alpha_trans))
cs.plot_cell_polygons(sham, "celltype"; 
    cell_colors=anno2, x_lims=(18300,19800), y_lims=(10700,14000), 
    width=300,height=600,stroke_color="gray80")

## inner medulla
anno2 = Dict("Podo" => ("magenta1",alpha_trans), 
                "HealthyPT"=>("green3",alpha_trans),
                "InjPT"=>("#f92874",alpha_trans),
                "TAL"=>("lightslateblue",alpha_trans),
                "DCT"=>("blue",alpha_trans),
                "CD-PC"=>("turquoise1",alpha_trans),
                "CD-IC"=>("#924cfa",alpha_trans),
                "vEC"=>("firebrick",alpha_trans),
                "gEC"=>("dodgerblue",alpha_trans),
                "Fib"=>("#edff4d",0.5),
                "JGA"=>("sienna2",alpha_trans),
                "Immune"=>("darkgreen",alpha_trans),
                "Uro"=>("black",alpha_trans))
cs.plot_cell_polygons(sham, "celltype"; 
    cell_colors=anno2, x_lims=(15000,17000), y_lims=(16000,20000), 
    width=300,height=600,stroke_color="gray40")

## papilla
anno2 = Dict("Podo" => ("magenta1",alpha_trans), 
                "HealthyPT"=>("green3",alpha_trans),
                "InjPT"=>("#f92874",alpha_trans),
                "TAL"=>("lightslateblue",alpha_trans),
                "DCT"=>("blue",alpha_trans),
                "CD-PC"=>("turquoise1",alpha_trans),
                "CD-IC"=>("#924cfa",alpha_trans),
                "vEC"=>("firebrick",alpha_trans),
                "gEC"=>("dodgerblue",alpha_trans),
                "Fib"=>("#edff4d",0.5),
                "JGA"=>("sienna2",alpha_trans),
                "Immune"=>("darkgreen",alpha_trans),
                "Uro"=>("black",alpha_trans));
cs.plot_cell_polygons(sham, "celltype"; 
    cell_colors=anno2, x_lims=(14300,16300), y_lims=(22900,24900), 
    width=400,height=400,stroke_color="gray40")




