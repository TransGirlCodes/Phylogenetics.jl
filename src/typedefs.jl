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








# Definitions for the clade variables.
# Type containing user defined properties, and the outer constructors.
type Property{T}
	attrRef::ASCIIString
	attrUnit::ASCIIString
	attrDatatype::DataType
	attrAppliesto::ASCIIString
	attrIdRef::ASCIIString
	value::T
end
# Properties outer constructor.
function Property{T}(v::T; ref::ASCIIString="", unit::ASCIIString="", appliesto::ASCIIString="", idref::ASCIIString="")
	typeval = typeof(v)
	Property(ref,unit,typeval,appliesto,idref,v)
end
# Type containing the confidence measure for a given clade feature.
type Confidence
	attrType::ASCIIString 
	value::Float64
end
# Type containing the uri information used in other types that clades contain.
type Uri
	attrDesc::ASCIIString
	attrType::ASCIIString
	string::ASCIIString
	Uri(a::ASCIIString,b::ASCIIString,c::ASCIIString) = new(a,b,c)
end
# Type containing the taxonomyID information for a clade taxonomy element.
type TaxonomyID
	provider::ASCIIString
	id::ASCIIString
	TaxonomyID(a::ASCIIString,b::ASCIIString) = new(a,b)
end
# Type containing taxonomy information for a Clade element.
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
	CladeTaxonomy(idsource::ASCIIString,id::TaxonomyID,code::ASCIIString,auth::ASCIIString,sname::ASCIIString,cname::Array{String},syn::Array{String},rank::ASCIIString,uri::Uri) = new(idsource,id,code,auth,sname,cname,syn,rank,uri)
end
# Type containing accession ifnroamtion for a sequence.
type SeqAccession
	source::ASCIIString
	accession::ASCIIString
	SeqAccession(source::ASCIIString,acc::ASCIIString) = new(source,acc)
end
# Type containing the molecular sequence.
type MolSeq
	aligned::Bool
	sequence::ASCIIString
	MolSeq(align::Bool, seq::ASCIIString) = new(align, seq)
end
# Type containing the information for a given annotation.
type SeqAnnotation
	attrRef::ASCIIString
	attrSource::ASCIIString
	attrEvidence::ASCIIString
	attrType::ASCIIString
	description::ASCIIString
	confidence::Confidence
	uri::Uri
	properties::Array{Property}
	SeqAnnotation(aref::ASCIIString, asource::ASCIIString, aevidence::ASCIIString, atype::ASCIIString, desc::ASCIIString, con::Confidence, prop::Array{Property}) = new(aref, asource, aevidence, atype, desc, con, prop)
end
# Type containing info for a protein domain.
type ProteinDomain
	attrFrom::Int64
	attrTo::Int64
	attrConfidence::Float64
	attrId::ASCIIString
	token::ASCIIString
	ProteinDomain(afrom::Int64, ato::Int64, aconf::Float64, aid::ASCIIString, token::ASCIIString) = new(afrom, ato, aconf, aid, token)
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
	proteinDomains::Array{ProteinDomain}
	CladeSequence(atype::ASCIIString, acc::SeqAccession, name::ASCIIString, symb::ASCIIString, mol::MolSeq, uri::Uri, ann::Array{SeqAnnotation}, dom::Array{ProteinDomain}) = new(atype, acc, name, symb, mol, uri, ann, dom)
end




















# Type containing the event information for any given clade.
type CladeEvents
	eventType::ASCIIString
	duplications::Int32
	speciations::Int32
	losses::Int32
	confidence::Array{Confidence}
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
	CladeColour(r::Float64, g::Float64, b::Float64) = new(r, g, b)
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

	