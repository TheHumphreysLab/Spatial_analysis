import CellScopes as cs

cellorder=["Podo","PTS12","PTS3","TAL","DCT","CD-PC","CD-IC","Uro"];
c_map=["magenta1","green3","olivedrab1","lightslateblue","blue","turquoise1","#924cfa","black"];
color1 = Dict(cellorder .=> c_map);
cs.sp_dim_plot(visium_male, :celltype; x_col="new_x", y_col="new_y", canvas_size=(800, 600), 
    do_label=false, marker_size=10, y_lims=(-5000, 1100), anno_color=color1, cell_order=cellorder)

genes=["Podxl","Slc5a2","Inmt","Csf1r","Havcr1","Umod","Slc12a3","Aqp2","Slc26a4","Krt19","Col1a1","Ren1","Cd74","Flt1","Eln"];
order_by=["Podo","PTS1","PTS2","PTS3","InjPT","TAL","DCT","CD-PC","CD-IC","Uro","Fib","JGA","Immune", "gEC","vEC"]
p3=cs.sp_dot_plot(sham,genes,:celltype2; 
    width=200,cell_order=order_by, height=250, expr_cutoff=1, perc_cutoff=20)

cellorder=["Podo","PTS12","PTS3","TAL","DCT","CD-PC","CD-IC","Uro"];
p4=sp_dot_plot(visium_sham,genes,:celltype; 
    width=200,cell_order=cellorder, height=250, expr_cutoff=1, perc_cutoff=20)
