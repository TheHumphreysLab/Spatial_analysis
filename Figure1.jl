using Pkg
using JLD2
import Baysor as B

pre_data=load("day2_baysor.jld2");
polygons= pre_data["polygons"];
df_spatial= pre_data["df_spatial"];

anno2 = Dict("1" => "#2e4ff2", "2"=>"yellow","3"=>"darkorange","4"=>"#03521f","5"=>"#b4b3b5","6"=>"#2cf51d",
            "7"=>"olivedrab1","8"=>"#8c49b8","9"=>"fuchsia","10"=>"#a100c9","11"=>"red","12"=>"cyan")

fig1b=B.plot_molecules(df_spatial, polygons, color=:cluster ,annotation=:cluster,
    markersize=1, size=(3000,3800),poly_strokewidth=0.5, ann_colors=anno2
)
MK.save("day2_plot_seg_whole.png", fig1)
fig1b

fig1c1=B.plot_molecules(df_spatial, polygons, color=:cluster ,annotation=:cluster,
    markersize=2, poly_strokewidth=1, ann_colors=anno2, size=(900,800),
    xlims=(21000,25000), ylims=(3600,6600)
);
MK.save("day2_plot_seg_cortex.png", fig1c1)
fig1c1

fig1c2=B.plot_molecules(df_spatial, polygons, color=:cluster ,annotation=:cluster,
    markersize=2, poly_strokewidth=1, ann_colors=anno2, size=(900,800),
    xlims=(16000,20500), ylims=(10000,13000)
);
MK.save("day2_plot_seg_om.png", fig1c2)
fig1c2

fig1c3=B.plot_molecules(df_spatial, polygons, color=:cluster ,annotation=:cluster,
    markersize=2, poly_strokewidth=1, ann_colors=anno2, size=(900,800),
    xlims=(15300,19300), ylims=(15200,18200)
);
MK.save("day2_plot_seg_im.png", fig1c3)
fig1c3

fig1c4=B.plot_molecules(df_spatial, polygons, color=:cluster ,annotation=:cluster,
    markersize=2, poly_strokewidth=1, ann_colors=anno2, size=(900,800),
    xlims=(13300,16300), ylims=(23000,26000)
);
MK.save("day2_plot_seg_p.png", fig1c4)
fig1c4

fig1c5=B.plot_molecules(df_spatial, polygons, color=:cluster ,annotation=:cluster,
    markersize=2, poly_strokewidth=1, ann_colors=anno2, size=(300,1000),
    xlims=(22000,24500), ylims=(11000,20000)
);
MK.save("day2_plot_seg_vasc.png", fig1c5)
fig1c5

