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
	attrIdSource::ASCIIString
	id::TaxonomyID
	code::ASCIIString
	authority::ASCIIString
	scientificName::ASCIIString
	commonName::Array{String}
	synonym::Array{String}
	rank::ASCIIString
	uri::Uri
end
type TaxonomyID
	provider::ASCIIString
	id::ASCIIString
end

# Type containing the uri information used in other types that clades contain.
type Uri
	attrDesc::ASCIIString
	attrType::ASCIIString
	string::ASCIIString
end

# Type containing the event information for any given clade.
type CladeEvents
	eventType::ASCIIString
	duplications::Int32
	speciations::Int32
	losses::Int32
	confidence::Confidence
end

# Type containing the confidence measure for a given clade feature.
type Confidence
	attrType::ASCIIString # The type of confidence value e.g. probability/bootstrap.
	value::Float64 # The actual value as a floating point number.
end

# Type and subtypes containing the sequence information for clades.
type CladeSequence
	attrType::ASCIIString
	accession::SeqAccession
	name::ASCIIString
	symbol::ASCIIString
	molSeq::MolSeq
	uri::Uri
	annotations::Array{SeqAnnotation}
	domainArchitecture::Array{}
end
type MolSeq
	aligned::Bool
	sequence::ASCIIString
end
type SeqAccession
	source::ASCIIString
	accession::ASCIIString
end
type SeqAnnotation
	attrRef::ASCIIString
	attrSource::ASCIIString
	attrEvidence::ASCIIString
	attrType::ASCIIString
	description::ASCIIString
	confidence::Confidence
	uri::Uri
	properties::Array{Property}
end

# Type containing user defined properties.
type Property
	attrRef::ASCIIString
	attrUnit::ASCIIString
	attrDatatype::ASCIIString
	attrAppliesto::ASCIIString
	attrIdRef::ASCIIString
	value
end



type BinaryCharacters
	attrType::ASCIIString
	attrGainedCount::Int64
	attrLostCount::Int64
	attrPresentCount::Int64
	attrAbsentCount::Int64
	gained::BinaryChar
	lost::BinaryChar
	present::BinaryChar
	absent::BinaryChar
end


type CladeColour
	red::Float64
	green::Float64
	blue::Float64
end




type PhyXClade
	attrIdSource::ASCIIString
	name::ASCIIString
	branchLength::Float64
	confidences::Array{Confidence}
	width::Float64
	branchColour::CladeColour
	nodeID::ASCIIString
	taxonomy::Array{CladeTaxonomy}
	sequence::PhyXSequence
	events::CladeEvents







end

	