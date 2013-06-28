module Phylo
	## Dependencies
	using Base.Intrinsics
	importall Base

	## Exported Methods and Types
	export Phylogeny,
		Phylo,
		Clado,
		ReducedTopology,
		TreeRead,
		treeWrite,
		getRoot,
		getKids

	## Load Package Files
	include(Pkg.dir("Phylo", "src", "typedefs.jl"))
	include(Pkg.dir("Phylo", "src", "treeio.jl"))

end

