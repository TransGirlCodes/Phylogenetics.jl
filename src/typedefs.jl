# Abstract type definition for Phylogenetic trees.
abstract Phylogeny

# The type definition for a Phylogenetic Tree with branch lengths.
type Phylo <: Phylogeny
	name::ASCIIString
	edge::Array{Int}
	Nnode::Int
	tipLabel::Array{ASCIIString}
	edgeLength::Array{Float64}
	nodeLabel::Array{String}
	rootEdge::Float64
	Phylo(name::ASCIIString, edge::Array{Int}, Nnode::Int, tipLabel::Array{ASCIIString}, edgeLength::Array{Float64}, nodeLabel::Array{String}, rootEdge::Float64) = new(name, edge, Nnode, tipLabel, edgeLength, nodeLabel, rootEdge)
end

# The type definition for a Phylogenetic Tree without branch lengths. 
type Clado <: Phylogeny
	name::ASCIIString
	edge::Array{Int}
	tipLabel::Array{ASCIIString}
	Nnode::Int
	nodeLabel::Array{String}
	Clado(name::ASCIIString, edge::Array{Int}, tipLabel::Array{ASCIIString}, Nnode::Int, nodeLabel::Array{String}) = new(name, edge, tipLabel, Nnode, nodeLabel)
end


# Type definition for a small simple representation of a tree. 
type ReducedTopology <: Phylogeny
	name::ASCIIString
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