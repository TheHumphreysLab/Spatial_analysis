using Pkg, StatsBase, CSV, DataFrames, VegaLite, Colors, ColorSchemes, JLD2
import CairoMakie as MK
import CellScopes as cs
using Statistics

# cell annotation for sham
sham=load("new_sham_SpaData_tg2.jld2")
sham=sham["sham"];
genes=["Podxl","Slc5a2","Inmt","Havcr1","Umod","Egf","Slc12a1","Slc12a3","Aqp2","Slc26a4","Flt1",
                      "Col1a1","Ren1","Acta2","Cd74","Krt19","Eln","Ehd3"];
cs.sp_dot_plot(sham,genes,:cluster; 
    width=300, height=300,  expr_cutoff=0.2)

cell1=collect(1:20)
cell2=["InjPT","TAL","vEC","HealthyPT","HealthyPT","Immune","TAL",
        "CD-PC","HealthyPT","gEC","Podo","TAL","TAL","Uro","TAL",
        "Fib","HealthyPT","JGA","DCT","CD-IC"]
sham.cells=cs.mapvalues(sham.cells, :cluster, :celltype, cell1, cell2);
cs.save(sham; filename="new_sham_SpaData_tg.jld2")

# cell annotation for hour4
hour4=load("new_hour4_SpaData_tg.jld2")
hour4=hour4["sp"];
cell1=string.(collect(1:20))
cell2=["HealthyPT","Podo","TAL","InjPT","InjPT","InjPT","HealthyPT",
        "gEC","HealthyPT","CD-IC","CD-PC","vEC","HealthyPT","Uro","TAL",
        "TAL","InjPT","Immune","Fib","DCT"]
hour4.cells=cs.mapvalues(hour4.cells, :cluster, :celltype, cell1, cell2);
cs.save(hour4; filename="new_hour4_SpaData_tg.jld2")

# cell annotation for hour12
hour12=load("new_hour12_SpaData_tg.jld2")
hour12=hour12["sp"];
cell1=string.(collect(1:20))
cell2=["TAL","JGA","TAL","TAL","InjPT","HealthyPT","InjPT",
        "gEC","TAL","HealthyPT","CD-IC","Podo","CD-PC","Fib","InjPT",
        "Immune","DCT","InjPT","vEC","InjPT"]
hour12.cells=cs.mapvalues(hour12.cells, :cluster, :celltype, cell1, cell2);
cs.save(hour12; filename="new_hour12_SpaData_tg.jld2")

# cell annotation for day2
day2=load("new_day2_SpaData_tg.jld2")
day2=day2["sp"];
cluster=collect(1:20)
celltypes=["Uro","DCT","InjPT","TAL","gEC","InjPT","Immune","JGA","Fib","InjPT","InjPT","vEC","HealthyPT",
    "Podo","InjPT","Fib","CD-IC","CD-PC","Fib","InjPT"]
day2.cells = cs.mapvalues(day2.cells, :cluster, :celltype, cluster, celltypes);
cs.save(day2; filename="new_day2_SpaData_tg.jld2")

# cell annotation for week6
week6=load("new_week6_SpaData_tg.jld2")
week6=week6["sp"];
cell1=collect(1:20)
cell2=["Uro","CD-IC","DCT","HealthyPT","Immune","CD-PC","CD-PC",
        "HealthyPT","Immune","InjPT","HealthyPT","Fib","JGA","vEC","TAL",
        "Immune","Immune","Podo","gEC","Fib"]
week6.cells=cs.mapvalues(week6.cells, :cluster, :celltype, cell1, cell2);
cs.save(week6; filename="new_week6_SpaData_tg.jld2")

# making the cell fraction df
sham_cell=cs.make_cell_proportion_df(sham.spmetaData.cell; nfeatures=0)
sham_cell.time.="Sham";
hour4_cell=cs.make_cell_proportion_df(hour4.spmetaData.cell; nfeatures=0)
hour4_cell.time.="Hour4";
hour12_cell=cs.make_cell_proportion_df(hour12.spmetaData.cell; nfeatures=0)
hour12_cell.time.="Hour12";
day2_cell=cs.make_cell_proportion_df(day2.cells; nfeatures=0)
day2_cell.time.="Day2";
week6_cell=cs.make_cell_proportion_df(week6.spmetaData.cell; nfeatures=0)
week6_cell.time.="Week6";
all_time=[sham_cell; hour4_cell; hour12_cell; day2_cell; week6_cell];
cell_order=["InjPT","HealthyPT","Podo","TAL","DCT","Fib","JGA",
            "vEC","gEC","CD-PC","CD-IC","Uro","Immune"];
cell_order2=["01-InjPT","02-HealthyPT","03-Podo","04-TAL","05-DCT","06-Fib","07-JGA",
            "08-vEC","09-gEC","10-CD-PC","11-CD-IC","12-Uro","13-Immune"];
cell_color=["#f92874","#5fca3b","#ea33f7","#8171f7","#0000cd","#edff4d","#df7f4f",
            "#eb463d","#3976cb" ,"#70f0fc","#924cfa" ,"#000000","#296218"];
all_time=cs.mapvalues(all_time, :celltype,:celltype2, cell_order, cell_order2);
CSV.write("cell_prop.csv", all_time)

#plot the cell fraction --figure 4A
all_time |> @vlplot(
    mark={:bar, opacity=0.7, width=50}, 
    x={:time, title="Time point", axis={labelFontSize=19,titleFontSize=19,grid=false},
         scale={domain=["Sham","Hour4","Hour12","Day2","Week6"]}}, 
    y={:fraction, title="Cell Fraction", axis={grid=false}},
    color={
        :celltype2,
        scale={
            domain=cell_order2,
            range=cell_color
        },
        legend={
            title="Cell type",titleFontSize=12,
            labelFontSize=12
        }
    },
    width=400,height=800)
time=["Sham","Hour4","Hour12","Day2","Week6"]
group_idx=["1-Sham","2-Hour4","3-Hour12","4-Day2","5-Week6"]
all_time=spd.mapvalues(all_time, :time, :index, time, group_idx)
p=all_time |> @vlplot()+@vlplot(mark={:area, opacity=0.6}, x={"index", axis={grid=false} }, y={:fraction, stack=:zero, axis={grid=false}}, color={"celltype2:n",  scale={
            domain=cell_order2,
            range=cell_color
        }},width=200,height=200) +
@vlplot(mark={:bar, width=1, opacity=1}, x={"index", title="Time point"}, y={:fraction, stack=:zero, title="Cell proportion"}, color={"celltype2:n",  scale={
            domain=cell_order2,
            range=cell_color
        }},
    width=200,height=200)
save("cell_proportion.pdf",p)

# plot gene across timepoints. repeat this for each timepoint and each gene.
cmap=ColorSchemes.ColorScheme([colorant"gray98",colorant"red", colorant"red4"])
day2=cs.generate_polygon_counts(day2);
fig_day2=cs.plot_gene_polygons(day2, "Col1a1",cmap; x_lims=(0, maximum(day2.cells.x)+1000), y_lims=(0, maximum(day2.cells.y)+1000),
    canvas_size=(6000,6000),stroke_width=0.5, stroke_color="gray60");
save("fig_day2_col1a1.pdf",fig_day2)



