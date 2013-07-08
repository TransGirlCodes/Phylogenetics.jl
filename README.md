Phylogenetics
=====
Version 0.0.1

**The Julia package for Analysis of Phylogeny and Evolution**

# Phylogenetics.jl
Phylogenetics.jl was started to mark the beginning of an explicit effort to bring means of working with phylogenetic data to Julia. 
It is still - like most of Julia - in development, however features are increasingly becoming useable.

# Installation 
```
Pkg.init() # Only the first time you install a Julia's Package

Pkg.add("Phylogenetics") # Install Phylogenetics.jl

using Phylogenetics # to use Phylogenetics
```

# Current Useable Features

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
		

Contributing
------------

**Fork and send a pull request or create a [GitHub issue](https://github.com/Ward9250/Phylo.jl/issues) for bug reports or feature requests.  Or if you created the feature you wanted put in a pull request!**