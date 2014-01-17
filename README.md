Phylogenetics
=====
Version 0.0.2

**The Julia package for Analysis of Phylogeny and Evolution**

# Installation
```
Pkg.init() # Only the first time you install a Julia's Package

Pkg.add("Phylogenetics") # Install Phylogenetics.jl

using Phylogenetics # to use Phylogenetics
```

# Features

* Types for Phylogenies
  * Abstract type "Phylogeny".
  * Type for Cladograms (trees without branch lengths) "Clado".
  * Type for phylogenetic trees with branch lengths "Phylo".
  * Type for reduced representation of trees as an array of indicies "ReducedTopology".
* Methods for reading in data from Newick format tree files.
  * Interprets information present in the newick tree and makes a Phylo or Clado type respectively. 
* Methods for writing out Phylo and Clado types as Newick format trees to file.
* Methods for manipulating Trees.
  * Find the root of the tree.
  * Display the children of every tree node.
  
# Features in Development / Requested

* Support for reading in phyloXML files (in branch phyxml).
  * A type for trees annotated with highly detailed information. 
  * Methods for converting the planned type with already existing types Clado, Phylo, and ReducedTopology.
		

Contributing
------------

**Fork and send a pull request or create a [GitHub issue](https://github.com/Ward9250/Phylo.jl/issues) for bug reports or feature requests.  Or if you created the feature you wanted put in a pull request! Request pulls to devel branch to keep master more stable and give time to find bugs before merge with master branch.**
