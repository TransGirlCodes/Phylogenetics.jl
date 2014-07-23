module Phylo
	## Dependencies - Base and Core are implictly imported!
  using LightXML
  import DataStructures: Queue, enqueue!, dequeue!, Stack

  ## Exported Methods and Types
	export PhyXExtension,
		PhyXElement,
		PhyXTree,

	## Load Package Files
	include(Pkg.dir("Phylogenetics", "src", "typedefs.jl"))
	include(Pkg.dir("Phylogenetics", "src", "treeio.jl"))
end
