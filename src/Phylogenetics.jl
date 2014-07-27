module Phylo

  ## Exported Methods and Types
  export PhyExtension,
    PhyNode,
    Phylogeny,
    getName,
    getBranchLength,
    isLeaf,
    hasChildren,
    parentIsSelf,
    hasParent,
    getChildren,
    getSiblings,
    getParent,
    isRoot,
    isNode,
    setName!,
    setBranchLength!,
    graft!,
    prune!

  ## Load Package Files
  include(Pkg.dir("Bio", "src", "phylo", "typedefs.jl"))
end
