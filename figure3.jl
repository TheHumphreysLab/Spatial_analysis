import CellScopes as cs
anno2 = Dict("Podo" => "magenta1", 
    "PTS1"=>"olivedrab1",
    "PTS2"=>"green3",
    "PTS3"=>"olive",
    "TAL"=>"lightslateblue",
    "DCT"=>"blue",
    "CD-PC"=>"turquoise1",
    "CD-IC"=>"#924cfa",
    "vEC"=>"firebrick",
    "gEC"=>"dodgerblue",
    "Fib"=>"#edff4d",
    "JGA"=>"sienna2",
    "Immune"=>"darkgreen",
    "Uro"=>"black");
cellorder=["Podo", "PTS1", "PTS2", "PTS3", "TAL", "DCT", "CD-PC", "CD-IC", "aEC", "gEC", "Fib", "JGA", "Immune", "Uro"];
cs.sp_dim_plot(male, :celltype; x_col = "x", y_col="y", canvas_size=(1200, 1100), 
    marker_size=2.6,stroke_width=0 , anno_color=anno2,label_size=40, 
    label_offset=(-0.3,0.5), do_legend=true, do_label=false, cell_order=cellorder)
cs.dim_plot(female, :celltype; x_col = "x", y_col="y", canvas_size=(1200, 1160), 
    marker_size=2.4,stroke_width=0 , anno_color=anno2,label_size=20, 
    label_offset=(-0.3,0.5), do_legend=true, do_label=false, cell_order=cellorder)
cs.sp_feature_plot(male,["Inmt","Rnf24", "Csf1r","Scel","Msln"]; marker_size=4 )
cs.sp_feature_plot(female,["Inmt","Rnf24", "Csf1r","Scel","Msln"]; marker_size=4 )

## heatmap
genes = ["Nphs2","Podxl", "Slc5a2","Slc5a12", "Inmt","Slc7a12","Tubb4b", "Csf1r","Rnf24", "Umod", "Egf", "Slc12a3", "Calb1", "Aqp2",
           "Slc26a4","Foxi1","Eln" ,"Tagln","Flt1","Emcn","Ehd3","Pecam1","Col1a1","Col6a1","Ren1","Cd74","Cd14", "Krt19","Msln"]
cs.plot_heatmap(male, genes, :celltype; assay_use = "measured", scale=true,
            fig_height=400, fig_width=150, color_scheme="darkgold", reverse_color=true,cell_order=cellorder)
cs.plot_heatmap(female, genes, :celltype; assay_use = "measured", scale=true,
            fig_height=400, fig_width=150, color_scheme="darkgold", reverse_color=true,cell_order=cellorder)
cellorder=["Podo", "PTS1", "PTS2", "PTS3", "TAL", "DCT", "CD-PC", "CD-IC", "vEC", "gEC", "Fib", "JGA", "Immune", "Uro"]
genes = ["Slc22a29","Dio1","Hao2","Prlr","Cyp4a14","Cyp2d26","Abcc3","Akr1c18","Serpina1f","Slc22a19",
    "Acsm3","Cml5","Cyp2e1","Cyp4a12a","Ugt2b38","Acsm3","Slc22a30","Cyp2j13","Fmo5","Azgp1","Slc7a13","Timd2","Ugt2b37"]
p1=cs.compare_gene_imputation(female, male, genes, :celltype; assay_use = "predicted",imp_type="SpaGE", scale=true,sp1_name="Female", sp2_name="Male", 
            fig_height=400, fig_width=150,cell_order=cellorder)
save("sex_specific_gene.pdf", p1)

## plot imputed genes
cs.sp_feature_plot(male,["Abcc3","Akr1c18", "Cyp2e1","Cyp4a12a"]; marker_size=4, use_imputed=true, imp_type = "SpaGE" )
cs.sp_feature_plot(female,["Abcc3","Akr1c18", "Cyp2e1","Cyp4a12a"]; marker_size=4, use_imputed=true, imp_type = "SpaGE" )




