Phylogenetics
=====
Version 0.0.2

**The Julia package for Analysis of Phylogeny and Evolution**

[![Phylogenetics](http://pkg.julialang.org/badges/Phylogenetics_0.2.svg)](http://pkg.julialang.org/?pkg=Phylogenetics&ver=0.2)

# Features in Development / Requested

* Current TODO list for Phylogenetics.jl is situated at [https://github.com/BioJulia/Bio.jl/issues/13](https://github.com/BioJulia/Bio.jl/issues/13)
* **Phylogenetics.jl is now being developed as a component Bio.jl, this repo should be viewed as an upstream/experimental repo of functionality that will be merged into the Phylo module of Bio.jl, assuming the code and ideas are accepted on review as fitting with a Bio.jl and Julian way of doing things.**


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
		

Contributing
------------

**Fork and send a pull request or create a [GitHub issue](https://github.com/Ward9250/Phylogenetics.jl/issues) for bug reports or feature requests.  Or if you created the feature you wanted put in a pull request! Request pulls to devel branch to keep master more stable and give time to find bugs before merge with master branch.**

MIT Liscence (c) 2014. Ben J. Ward.
