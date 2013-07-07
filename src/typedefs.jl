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

# Definition of PhyXTree and it's elements.
type PhyXTree
	name::ASCIIString
	structure::Array{Int}
	tipNodes::Array{Int}
	internalNodes::Array{Int}
	userData::Array{Array}
end

# Definition for the clade elements.


# Definitions for the clade variables
# Type containing the taxonomy information for a clade element.
type CladeTaxonomy
	provider::ASCIIString
	taxonomyID::ASCIIString
	taxonomyCode::ASCIIString
	taxonomyAuthority::ASCIIString
	scientificName::ASCIIString
	commonName::Array{String}
	synonym::Array{String}
	rank::ASCIIString
	uriType::ASCIIString
	uriDesc::ASCIIString
	uri::ASCIIString
end





type CladeSequence
	database::ASCIIString
	accession::ASCIIString
	name::ASCIIString
	symbol::ASCIIString
	molSeq::MolSeq
	binarySeq::ASCIIString
end

type PhyXClade
	cladeID::ASCIIString
	taxonomy::PhyXTaxonomy
	sequence::PhyXSequence
	events::CladeEvents
end

type CladeEvents
	speciations::Int
end






type MolSeq
	aligned::Bool
	sequence::ASCIIString
end

	