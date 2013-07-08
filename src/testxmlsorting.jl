# Changing approach to PhyloXML parsing useing libexpat.


using LibExpat
filepath = "~/Desktop/apaf-1.xml"
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
	function process_uri(instr)
		uritype = instr[xpath"uri/@type"]
		uridesc = instr[xpath"uri/@desc"]
		uri = instr[xpath"uri/text()"]
		uritype = length(uritype) > 0 ? uritype[1] : ""
		uridesc = length(uridesc) > 0 ? uridesc[1] : ""
		uri = length(uri) > 0 ? uri[1] : ""
		outUri = Uri(uridesc, uritype, uri)
		return outUri
	end
	function process_taxonomy(inclade)
		# Extract the values out of the elements.
		taxonomies = inclade[xpath"taxonomy"]
		outTaxonomies = Array(CladeTaxonomy, length(taxonomies))
		for i in 1:length(taxonomies)
			taxo = taxonomies[i]
			attrIdSource = taxo[xpath"@id_source"]
			idProvider = taxo[xpath"id/@provider"]
			taxonomyId = taxo[xpath"id/text()"]
			taxonomyCode = taxo[xpath"code/text()"]
			scientificName = taxo[xpath"scientific_name/text()"]
			taxonomyAuthority = taxo[xpath"authority/text()"]
			commonName = taxo[xpath"common_name/text()"]
			synonym = taxo[xpath"synonym/text()"]
			rank = taxo[xpath"rank/text()"]
			outUri = process_uri(taxo)
			# Sort out the values retrieved.
			attrIdSource = length(attrIdSource) > 0 ? attrIdSource[1] : ""
			scientificName = length(scientificName) > 0 ? scientificName[1] : ""
			idProvider = length(idProvider) > 0 ? idProvider[1] : ""
			taxonomyId = length(taxonomyId) > 0 ? taxonomyId[1] : ""
			taxonomyCode = length(taxonomyCode) > 0 ? taxonomyCode[1] : ""
			taxonomyAuthority = length(taxonomyAuthority) > 0 ? taxonomyAuthority[1] : ""
			rank = length(rank) > 0 ? rank[1] : ""
			# Make the taxonomy type.
			outid = TaxonomyID(idProvider, taxonomyId)
			outTaxonomies[i] = CladeTaxonomy(attrIdSource, outid, taxonomyCode, taxonomyAuthority, scientificName, commonName, synonym, rank, outUri)
		end
		return outTaxonomies
	end




	function process_confidence(instr)
		confidence = instr[xpath"number()"]
		confidenceType = instr[xpath"@type"]
		confidenceType = length(confidenceType) > 0 ? confidenceType[1] : ""
		outConfidence = Confidence(confidenceType, confidence)
		return outConfidence
	end
	
	function process_events(inclade)
		# Get the values out of the elements.
		eventType = inclade[xpath"events/type/text()"]
		duplications = int32(inclade[xpath"events/duplications/text()"])
		speciations = int32(inclade[xpath"events/speciations/text()"])
		losses = int32(inclade[xpath"events/losses/text()"])
		outConfidence = process_confidence(inclade[xpath"events"])[1]
		# Sort out the values retrieved.
		eventType = length(eventType) > 0 ? eventType[1] : ""
		duplications = length(duplications) > 0 ? duplications[1] : int32(0)
		speciations = length(speciations) > 0 ? speciations[1] : int32(0)
		losses = length(losses) > 0 ? losses[1] : int32(0)
		# Make the events type.
		outEvents = CladeEvents(eventType, duplications, speciations, losses, outConfidence)
		return outEvents
	end
	
	function process_properties(instr)
		properties = instr[xpath"property"]
		propertiesOut = Array(Property, length(properties))
		for i in 1:length(properties)
			property = properties[i]
			attrref = property[xpath"@ref"]
			attrunit = property[xpath"@unit"]
			attrdatatype = property[xpath"@datatype"]
			attrappliesto = property[xpath"@applies_to"]
			attridref = property[xpath"@id_ref"]
			value = property[xpath"text()"]
			attrref = length(attrref) > 0 ? attrref[1] : ""
			attrunit = length(attrunit) > 0 ? attrunit[1] : ""
			attrdatatype = length(attrdatatype) > 0 ? attrdatatype[1] : ""
			attrappliesto = length(attrappliesto) > 0 ? attrappliesto[1] : ""
			attridref = length(attridref) > 0 ? attridref[1] : ""
			propertiesOut[i] = Properties(attrref, attrunit, attrdatatype, attrappliesto, attridref, value)
		end
		return propertiesOut
	end

	function process_annotations(instr)
		annotationsOut = Array(SeqAnnotation, length(instr))
		for i in 1:length(annotationsOut)
			annotation = instr[i]
			attrref =  annotation[xpath"@ref"]
			attrsource = annotation[xpath"@source"]
			attrevidence = annotation[xpath"@evidence"]
			attrtype = annotation[xpath"@type"]
			desc = annotation[xpath"desc/text()"]
			attrref = length(attrref) > 0 ? attrref[1] : ""
			attrsource = length(attrsource) > 0 ? attrsource[1] : ""
			attrevidence = length(attrevidence) > 0 ? attrevidence[1] : ""
			attrtype = length(attrtype) > 0 ? attrtype[1] : ""
			desc = length(desc) > 0 ? desc[1] : ""
			confidence = process_confidence(annotation)
			uri = process_uri(annotation)
			properties = process_properties(annotation)
			annotationsOut[i] = SeqAnnotation(attrref, attrsource, attrevidence, attrtype, desc, confidence, uri, properties)
		end
		return annotationsOut
	end
	function process_sequence(inclade)
		sequences = inclade[xpath"sequence"]
		outSequences = Array(CladeSequence, length(sequences))
		for i in 1:length(sequences)
			seq = sequences[i]
			# Extract the values out of the elements.
			attrType = seq[xpath"sequence/@type"]
			seqAccessionSource = seq[xpath"accession/@source"]
			seqAccession = seq[xpath"accession/text()"]
			seqName = seq[xpath"name/text()"]
			molSeq = seq[xpath"mol_seq/text()"]
			molSeqAligned = seq[xpath"mol_seq/@is_aligned"]
			symbol = seq[xpath"symbol/text()"]
			# Sort out the values returned.
			attrType = length(attrType) > 0 ? attrType[1] : ""
			seqAccessionSource = length(seqAccessionSource) > 0 ? seqAccessionSource[1] : ""
			seqAccession = length(seqAccession) > 0 ? seqAccession[1] : ""
			seqName = length(seqName) > 0 ? seqName[1] : ""
			molSeq = length(molSeq) > 0 ? molSeq[1] : ""
			molSeqAligned = length(molSeqAligned) > 0 ? molSeq[1] : false
			symbol = length(symbol) > 0 ? symbol[1] : ""
			# Make the sequence type.
			annotations = process_annotations(seq[xpath"annotation"])
			uriOut = process_uri(seq)
			outAccession = SeqAccession(seqAccessionSource, seqAccession)
			molecularSequence = MolSeq(molSeqAligned, molSeq)
			domainarch = process_domainArchitecture(seq)
			outSequences[i] = CladeSequence(attrType, outAccession, seqName, symbol, molecularSequence, uriOut, annotations)
		end
		return outSequences
	end
	
	function process_binary(inclade)
		binary = inclade[xpath"binary_characters"]

	end

	function process_domainArchitecture(instr)
		domains = instr[xpath"domain_architecture/domain"]
		attrlength = domains[xpath"domain_architecture/@length"]
		attrlength = length(attrlength) > 0 ? attrlength[1] : 0
		outDomains = Array(Domain, length(domains))
		for i in 1:length(domains)
			dom = domains[i]
			attrfrom = dom[xpath"@from"]
			attrto = dom[xpath"@to"]
			attrconfidence = dom[xpath"@confidence"]
			attrid = dom[xpath"@id"]
			token = dom[xpath"text()"]
			attrfrom = length(attrfrom) > 0 ? attrfrom[1] : -1
			attrto = length(attrto) > 0 ? attrto[1] : -1
			attrconfidence = length(attrconfidence) > 0 ? attrconfidence[1] : -1
			attrid = length(attrid) > 0 ? attrid[1] : -1
			attrfrom = length(attrfrom) > 0 ? attrfrom[1] : -1
			outDomains[i] = Domain(attrfrom, attrto, attrconfidence, attrid, token)
		end


	end

	










# If no number is present in the xpath you've set it will return NaN or not a number.


	rooted = get(tree.attr, "rooted", "") == "true" ? true : false
	rerootable = get(tree.attr, "rerootable", "") == "true" ? true : false
	nClades = length(find(tree, "//clade"))
	allClades = find(tree, "//clade")
	structure = zeros(Int, nClades)
	outClades = Array(PhyXClade, nClades)	

	for i in 1:length(allClades)
		currentClade = allClades[i]
		children = find(currentClade, "clade")
		if length(children) == 0
			cladeType = 1 # Clade Type 1 is a tip.
		elseif length(children) > 0
			cladeType = 2 # Clade Type 2 is an internal node.
			childnodes = findin(allClades, children)
			structure[childnodes] = i
		end
		name = currentClade[xpath"name/text()"]
		name = length(name) > 0 ? name[1] : ""
		branchLength = currentClade[xpath"number(@branch_length)"]
		if isnan(branchLength) && !isnan(currentClade[xpath"number(branch_length)"])
			branchLength = currentClade[xpath"number(branch_length)"]
		end
		confs = currentClade[xpath"confidence"]
		confidences = Array(Confidence, length(confs))
		for i in 1:length(confs)
			confidences[i] = process_confidence(confs[i])
		end

		
		# All these things appear to work.
		width = currentClade[xpath"number(width)"]
		colr = currentClade[xpath"number(color/red)"]
		colg = currentClade[xpath"number(color/green)"]
		colb = currentClade[xpath"number(color/blue)"]
		branchcol = CladeColour(colr, colg, colb)
		nodeid = currentClade[xpath"node_id/text()"]
		nodeid = length(nodeid) > 0 ? nodeid[1] : ""
		taxonomies = process_taxonomy(currentClade)



		sequence = process_sequence(currentClade)
		


		events = process_events(currentClade)
		binarychars = process_binary(currentClade)


		
		
		
		
		



	end





	edge = hcat(structure, [1:nClades])
	ind = [findin(structure, i) for i in [1,2,3,4,5]] # Figure which nodes have kids.
	tips = [1:nClades][[i == [] for i in ind]] # Tip nodes should have no kids.
	nodes = [1:nClades][[i != [] for i in ind]] # Internal nodes should have kids.
	newNodes = [length(tips) + i for i in 1:length(nodes)]
end