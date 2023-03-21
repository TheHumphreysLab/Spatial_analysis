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

using RCall
@rlibrary ggplot2
@rlibrary scatterpie
vs_data = deepcopy(visium_male.spmetaData.cell)
sp_data = deepcopy(cartana_male.spmetaData.cell)
ggplt=rimport("ggplot2")
cellorder1=["Podo","PTS1","PTS2","PTS3","TAL","DCT","CD-PC","CD-IC","vEC","gEC", "Fib","JGA","Immune", "Uro"]
colors1=["magenta1","olivedrab1","green","darkolivegreen4","lightslateblue","blue","turquoise1","purple","firebrick",
"dodgerblue","#edff4d","sienna2","darkgreen","black"];
vs_r = visium_unit_radius(41)
p = ggplot(sp_data, aes(:new_x, :new_y))+geom_point(size=6, color="gray85", alpha=0.5)+ 
ggplt.xlim(0.57, 1) + ggplt.ylim(0.18, 0.5)+
scale_color_manual(values=colors1)+
ggplt.theme(var"panel.background" = element_rect(fill = "white", colour = "black"),
          var"axis.ticks" = element_blank(),
          var"axis.line" = element_blank(),
          var"axis.text" = element_blank(),
          var"axis.title" =element_blank(),
          var"strip.text" = element_text(size=12),
          var"legend.title" = element_blank(),
          var"legend.text"=element_text(size=40),
          var"legend.key.size" = unit(2, "cm"),
          var"legend.position" ="right")
p1=p + geom_scatterpie(aes(x=:new_x, y=:new_y, group=:cell, r=vs_r),
                    data=vs_data, cols=cellorder1, alpha=1) +
scale_fill_manual(values=colors1)

sp_data=cs.reorder(sp_data, :celltype, cellorder1)
colors2=["purple","turquoise1","blue", "#edff4d","dodgerblue","darkgreen","sienna2","magenta1",
    "olivedrab1","green","darkolivegreen4", "lightslateblue", "black","firebrick"]
ggplot(sp_data, aes(x=:new_x, y=:new_y, color=:celltype))+geom_point(size=6, alpha=1)+ 
ggplt.xlim(0.57, 1) + ggplt.ylim(0.18, 0.5)+
scale_color_manual(values=colors2)+
ggplt.theme(var"panel.background" = element_rect(fill = "white", colour = "black"),
          var"axis.ticks" = element_blank(),
          var"axis.line" = element_blank(),
          var"axis.text" = element_blank(),
          var"axis.title" =element_blank(),
          var"strip.text" = element_text(size=12),
          var"legend.title" = element_blank(),
          var"legend.text"=element_text(size=40),
        var"legend.key" = element_rect(fill = "transparent"),
          var"legend.key.size" = unit(2, "cm"),
          var"legend.position" ="right")+
 guides(color = guide_legend(var"override.aes" = R"list"(size = 18)))
