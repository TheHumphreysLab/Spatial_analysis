using Pkg, StatsBase, CSV, DataFrames, VegaLite, Colors, ColorSchemes, JLD2
import CairoMakie as MK
import SpaData as spd
using PyCall
using Downloads
ENV["PYTHON"]="/bin/miniconda3/bin/python"
Pkg.build("PyCall")
sc = pyimport("scanpy")

Downloads.download("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139107&format=file&file=GSE139107%5FMouseIRI%2Emetadata%2Etxt%2Egz", "IRI_metadata.txt.tar.gz")
sc_meta=CSV.File("IRI_metadata.txt.tar.gz"; header=2) |> DataFrame
Downloads.download("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139107&format=file&file=GSE139107%5FMouseIRI%5Fcontrol%2Edge%2Etxt%2Egz", "sham.txt.tar.gz")
sc_count=CSV.File("sham.txt.tar.gz"; header=2) |> DataFrame
