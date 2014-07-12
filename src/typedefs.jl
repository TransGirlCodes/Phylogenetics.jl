## Recurvive extendible type for representation of phylogenetic trees in in Julia.

type PhyXExtension{T}
  value::T
end

type PhyXElement
  Label::String
  Root::Bool
  Tip::Bool
  Extensions::Array{PhyXExtension}
  Parent::PhyXElement
  Children::Array{PhyXElement}

  PhyXElement(label::String, root::Bool, tip::Bool, ext::Array{PhyXExtension, 1}) = new(label, root, tip, ext)
end

function PhyXElement(label::String, root::Bool, tip::Bool, ext::Array{PhyXExtension, 1}, parent::PhyXElement)
  x = PhyXElement(label, root, tip, ext)
  x.Parent = parent
  return x
end

type PhyXTree
  Name::String
  Root::PhyXElement
  Rooted::Bool
  Rerootable::Bool
end

# Abstract type definition for Phylogenetic trees.
abstract Phylogeny

# The type definition for a Phylogenetic Tree with branch lengths.
immutable Phylo <: Phylogeny
	name::String
	edge::Array{Int,2}
	Nnode::Int
	tipLabel::Array{String}
	edgeLength::Array{Float64}
	nodeLabel::Array{String}
	rootEdge::Float64
	Phylo(name, edge,
              Nnode, tipLabel,
              edgeLength, nodeLabel,
              rootEdge) = new(name, edge,
                              Nnode, tipLabel,
                              edgeLength, nodeLabel, rootEdge)
end

# The type definition for a Phylogenetic Tree without branch lengths.
immutable Clado <: Phylogeny
	name::String
	edge::Array{Int,2}
	Nnode::Int
	tipLabel::Array{String}
	nodeLabel::Array{String}
	Clado(name,
              edge,
              Nnode,
              tipLabel,
              nodeLabel) = new(name, edge, Nnode, tipLabel, nodeLabel)
end

# Type definition for a small simple representation of a tree.
immutable ReducedTopology <: Phylogeny
	name::String
	indiesArray::Array{Int}

	function ReducedTopology(phy::Phylogeny)
		children = phy.edge[1:size(phy.edge,1), 2]
		parents = phy.edge[1:size(phy.edge,1), 1]
		IArray = [length(findin(children, i)) != 0 ? parents[findin(children, i)[1]] : -1 for i in 1:max(phy.edge)]
		new(phy.name, IArray)
	end
end

function ReducedTopology(phy::Array{Phylogeny})
	outarray = Array(ReducedTopology, length(phy))
	for i in 1:length(phy)
		outarray[i] = ReducedTopology(phy[i])
	end
	return outarray
end

# equality, etc.

function =={T<:Phylogeny}(p1::T,p2::T)
    for field in names(T)
        if getfield(p1,field) != getfield(p2,field)
            return false
        end
    end
    return true
end

# Julia 0.2 compatibility
function isequal{T<:Phylogeny}(p1::T,p2::T)
    for field in names(T)
        if !isequal(getfield(p1,field),getfield(p2,field))
            return false
        end
    end
    return true
end

function hash{T<:Phylogeny}(p::T)
    h = 0
    for field in names(T)
        h = bitmix(hash(getfield(p,field)),h)
    end
    return h
end
