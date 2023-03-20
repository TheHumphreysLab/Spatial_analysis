import CellScopes as cs
sham = cs.load(filename="sham_final.jld2")
hour4 = cs.load(filename="hour4_final.jld2")
hour12 = cs.load(filename="hour12_final.jld2")
day2 = cs.load(filename="day2_final.jld2")
week6 = cs.load(filename="week6_final.jld2")

p1=cs.highlight_cells(sham, "InjPT", :celltype; cell_color="#f92874", stroke_color="gray80", marker_size=4);
p2=cs.highlight_cells(hour4, "InjPT", :celltype; cell_color="#f92874", stroke_color="gray80", marker_size=4);
p3=cs.highlight_cells(hour12, "InjPT", :celltype; cell_color="#f92874", stroke_color="gray80", marker_size=4);
p4=cs.highlight_cells(day2, "InjPT", :celltype; cell_color="#f92874", stroke_color="gray80", marker_size=4);
p5=cs.highlight_cells(week6, "InjPT", :celltype; cell_color="#f92874", stroke_color="gray80", marker_size=4);
save("injPT_p1.pdf", p1)
save("injPT_p2.pdf", p2)
save("injPT_p3.pdf", p3)
save("injPT_p4.pdf", p4)
save("injPT_p5.pdf", p5);

p1=cs.highlight_cells(sham, "Fib", :celltype; cell_color="gold1", stroke_color="gray90", marker_size=4);
p2=cs.highlight_cells(hour4, "Fib", :celltype; cell_color="gold1", stroke_color="gray90", marker_size=4);
p3=cs.highlight_cells(hour12, "Fib", :celltype; cell_color="gold1", stroke_color="gray90", marker_size=4);
p4=cs.highlight_cells(day2, "Fib", :celltype; cell_color="gold1", stroke_color="gray90", marker_size=4);
p5=cs.highlight_cells(week6, "Fib", :celltype; cell_color="gold1", stroke_color="gray90", marker_size=4);
save("fib_p1.pdf", p1)
save("fib_p2.pdf", p2)
save("fib_p3.pdf", p3)
save("fib_p4.pdf", p4)
save("fib_p5.pdf", p5);

p1=cs.highlight_cells(sham, "Immune", :celltype; cell_color="darkgreen", stroke_color="gray80", marker_size=4);
p2=cs.highlight_cells(hour4, "Immune", :celltype; cell_color="darkgreen", stroke_color="gray80", marker_size=4);
p3=cs.highlight_cells(hour12, "Immune", :celltype; cell_color="darkgreen", stroke_color="gray80", marker_size=4);
p4=cs.highlight_cells(day2, "Immune", :celltype; cell_color="darkgreen", stroke_color="gray80", marker_size=4);
p5=cs.highlight_cells(week6, "Immune", :celltype; cell_color="darkgreen", stroke_color="gray80", marker_size=4);
save("Immune_p1.pdf", p1)
save("Immune_p2.pdf", p2)
save("Immune_p3.pdf", p3)
save("Immune_p4.pdf", p4)
save("Immune_p5.pdf", p5);
