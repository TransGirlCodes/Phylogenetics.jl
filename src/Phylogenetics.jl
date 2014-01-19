module Phylogenetics
	## Dependencies
	using Base.Intrinsics
	importall Base

	## Exported Methods and Types
	export Phylogeny,
		Phylo,
		Clado,
		ReducedTopology,
		readtree,
		writetree,
		getroot,
		getkids,
		@tr_str

	## Load Package Files
	include(Pkg.dir("Phylogenetics", "src", "typedefs.jl"))
	include(Pkg.dir("Phylogenetics", "src", "treeio.jl"))

end
