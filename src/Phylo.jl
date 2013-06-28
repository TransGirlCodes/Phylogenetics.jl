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
	include(Pkg.dir("Jape", "src", "typedefs.jl"))
	include(Pkg.dir("Jape", "src", "treeio.jl"))

end

