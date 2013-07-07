# Changing approach to PhyloXML parsing useing libexpat.


using LibExpat
filepath = "~/Desktop/phyxml"
instream = open(expanduser(filepath))
instring = readall(instream)
close(instream)
xmltree = xp_parse(instring)
xmltree[xpath"phyloxml/phylogeny/clade"]



function phyloxbuild(inputstring)
	xmltree = xp_parse(inputstring)
	allInTrees = xmltree[xpath"phylogeny"]
	if length(allInTrees) == 1
		OutTree = buildphyx(allInTrees[1])
		return OutTree
	end
	if length(allInTrees) > 1
		allOutTrees = Array(Any, length(allInTrees))
		for tree in 1:length(alltrees)
			allOutTrees[tree] = buildphyx(allInTrees[tree])
		end
		return allOutTrees
	end
end



function buildphyx(tree)
	function process_taxonomy(inclade)
		# Extract the values out of the elements.
		idProvider = inclade[xpath"taxonomy/id/@provider"]
		taxonomyId = inclade[xpath"taxonomy/id/text()"]
		taxonomyCode = inclade[xpath"taxonomy/code/text()"]
		scientificName = inclade[xpath"taxonomy/scientific_name/text()"]
		taxonomyAuthority = inclade[xpath"taxonomy/authority/text()"]
		commonName = inclade[xpath"taxonomy/common_name/text()"]
		synonym = inclade[xpath"taxonomy/synonym/text()"]
		rank = inclade[xpath"taxonomy/rank/text()"]
		uritype = inclade[xpath"taxonomy/uri/@type"]
		uridesc = inclade[xpath"taxonomy/uri/@desc"]
		uri = inclade[xpath"taxonomy/uri/text()"]
		# Sort out the values retrieved.
		scientificName = length(scientificName) > 0 ? scientificName[1] : ""
		idProvider = length(idProvider) > 0 ? idProvider[1] : ""
		taxonomyId = length(taxonomyId) > 0 ? taxonomyId[1] : ""
		taxonomyCode = length(taxonomyCode) > 0 ? taxonomyCode[1] : ""
		taxonomyAuthority = length(taxonomyAuthority) > 0 ? taxonomyAuthority[1] : ""
		rank = length(rank) > 0 ? rank[1] : ""
		uritype = length(uritype) > 0 ? uritype[1] : ""
		uridesc = length(uridesc) > 0 ? uridesc[1] : ""
		uri = length(uri) > 0 ? uri[1] : ""
		# Make the taxonomy type.
		outTaxonomy = CladeTaxonomy(idProvider, taxonomyId, taxonomyCode, taxonomyAuthority, scientificName, commonName, rank, uritype, uridesc, uri)
		return outTaxonomy
	end



	function process_sequence(inclade)
		# Extract the values out of the elements.
		seqAccessionSource = inclade[xpath"sequence/accession/@source"]
		seqAccession = inclade[xpath"sequence/accession/text()"]
		seqName = inclade[xpath"sequence/name/text()"]
		molSeq = inclade[xpath"sequence/mol_seq/text()"]
		molSeqAligned = inclade[xpath"sequence/mol_seq/@is_aligned"]
		symbol = inclade[xpath"sequence/symbol/text()"]
		binSeq = 
		# Sort out the values retrieved.
		seqAccessionSource = length(seqAccessionSource) > 0 ? seqAccessionSource[1] : ""
		seqAccession = length(seqAccession) > 0 ? seqAccession[1] : ""
		seqName = length(seqName) > 0 ? seqName[1] : ""
		molSeq = length(molSeq) > 0 ? molSeq[1] : ""
		molSeqAligned = length(molSeqAligned) > 0 ? molSeq[1] : false
		symbol = length(symbol) > 0 ? symbol[1] : ""
		# Make the sequence type.
		molecularSequence = MolSeq(molSeqAligned, molSeq)
		outSequence = CladeSequence(seqAccessionSource, seqAccession, seqName, symbol, molecularSequence)
	end

	rooted = get(tree.attr, "rooted", "")
	rerootable = get(tree.attr, "rerootable", "")
	nClades = length(find(tree, "//clade"))
	allClades = find(tree, "//clade")
	structure = zeros(Int, nClades)
	for i in 1:length(allClades)
		currentClade = allClades[i]
		children = find(currentClade, "clade")
		if length(children) == 0
			cladeType = 1 # Clade Type 1 is a tip.
			taxonomy = process_taxonomy(currentClade)
			sequence = process_sequence(currentClade)


		elseif length(children) > 0
			cladeType = 2 # Clade Type 2 is an internal node.
			childnodes = findin(allClades, children)
			structure[childnodes] = i
			currentClade[xpath"number(/events/speciations)"]
			speciations = int32(currentClade[xpath"events/speciations/text()"])			
		end
		


	end





	edge = hcat(structure, [1:nClades])
	ind = [findin(structure, i) for i in [1,2,3,4,5]] # Figure which nodes have kids.
	tips = [1:nClades][[i == [] for i in ind]] # Tip nodes should have no kids.
	nodes = [1:nClades][[i != [] for i in ind]] # Internal nodes should have kids.
	newNodes = [length(tips) + i for i in 1:length(nodes)]
end