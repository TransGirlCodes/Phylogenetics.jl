



# Function method for writing a single Phylogeny type to file.
function treewrite(tree::Phylogeny,
                   file::ASCIIString = "output.nwk",
                   append::Bool = false,
                   treeNames::Bool = true)
	output = newick(tree, treeNames)
	if append == true
		outstream = open(file, "a")
	else
		outstream = open(file, "w")
	end
	println(outstream, output)
	close(outstream)
end


# Function method for writing multiple trees to file. In order to do this they must be as an array.
function treewrite(tree::Array{Phylogeny},
                   file::ASCIIString = "output.nwk",
                   append::Bool = false,
                   treeNames::Bool = false)
	outputarray = Array(ASCIIString, length(tree))
	for i in 1:length(tree)
		outputarray[i] = newick(tree[i], treeNames)
	end
	if apppend == true
		outstream = open(file, "a")
	else
		outstream = open(file, "w")
	end
	for i in outputarray
		println(outstream, i)
	end
	close(outstream)
end







# This function returns an array of arrays, showing which children
# nodes have.
function getkids(phy::Phylogeny)
	N = length(phy.tipLabel)
	kids = Array(Array{Int64}, N + phy.Nnode)
	for i in 1:N + phy.Nnode
		logic = bool([n == i for n in phy.edge[1:size(phy.edge,1), 1]])
		kids[i] = phy.edge[logic,2]
	end
	return kids
end


# Construct a newick string from a Cladogram
function newick(phy::Clado, name::Bool)
	function addInternal(i, k, STRING, N, nodelab, tiplab, ind)
		k, STRING = cp("(", k, STRING)
		desc = kids[i]
		for j in desc
			if j > N
				STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab)
			else
				k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			end
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
		k, STRING = cp(")", k, STRING)
		if i > N
			k, STRING = cp(nodelab[i - N], k, STRING)
		end
		return STRING, k
	end
	function addTerminal(i, k, STRING, tiplab, children)
		i = i[1]
		k, STRING = cp(tiplab[children[i]], k, STRING)
		return k, STRING
	end
	if name == true
		prefix = phy.name
	else
		prefix = ""
	end
	nodelab = [i != "" ? checklabels(i) : "" for i in phy.nodeLabel]
	tiplab = [i != "" ? checklabels(i) : i for i in phy.tipLabel]
	children = phy.edge[1:size(phy.edge,1), 2]
	parents = phy.edge[1:size(phy.edge,1), 1]
	N = length(phy.tipLabel)
	kids = getkids(phy)
	LS = (4 * N) + 5
	LS = LS + N # if there are nodelabels.
	ind = [findin(children, i) for i in 1:maximum(phy.edge)]
	STRING = ["" for i in 1:LS]
	k = 1
	k, STRING = cp(prefix, k, STRING)
	k, STRING = cp("(", k, STRING)
	desc = kids[getroot(children, parents)]
	for j in desc
		if j > N
			STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab, ind)
		else
			k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
	end
	k, STRING = cp(")", k, STRING)
	k, STRING = cp(nodelab[1], k, STRING)
	k, STRING = cp(";", k, STRING)
	outstring = ""
	for i in STRING
		outstring = "$outstring$i"
	end
	if name == true && prefix != ""
		namebit = chop(match(r"^[^\(]+\(", outstring).match)
		replace(outstring, namebit, "$namebit ")
	end
	return outstring
end


# Function that creates
function newick(phy::Phylo, name::Bool)
	function addInternal(i, k, STRING, N, nodelab, tiplab, ind)
		k, STRING = cp("(", k, STRING)
		desc = kids[i]
		for j in desc
			if j > N
				STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab)
			else
				k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			end
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
		k, STRING = cp(")", k, STRING)
		if i > N
			k, STRING = cp(nodelab[i - N], k, STRING)
		end
		k, STRING = cp(":", k, STRING)
		edgechar = phy.edgeLength[ind[i]][1]
		k, STRING = cp("$edgechar", k, STRING)
		return STRING, k
	end
	function addTerminal(i, k, STRING, tiplab, children)
		i = i[1]
		k, STRING = cp(tiplab[children[i]], k, STRING)
		k, STRING = cp(":", k, STRING)
		edgechar = phy.edgeLength[i]
		k, STRING = cp("$edgechar", k, STRING)
		return k, STRING
	end
	if name == true
		prefix = phy.name
	else
		prefix = ""
	end
	brl = phy.edgeLength
	nodelab = [i != "" ? checklabels(i) : "" for i in phy.nodeLabel]
	tiplab = [i != "" ? checklabels(i) : i for i in phy.tipLabel]
	children = phy.edge[1:size(phy.edge,1), 2]
	parents = phy.edge[1:size(phy.edge,1), 1]
	N = length(phy.tipLabel)
	kids = getkids(phy)
	LS = (4 * N) + 5
	LS = LS + N
	LS = LS +(4 * N)
	ind = [findin(children, i) for i in 1:maximum(phy.edge)]
	STRING = ["" for i in 1:LS]
	k = 1
	k, STRING = cp(prefix, k, STRING)
	k, STRING = cp("(", k, STRING)
	desc = kids[getroot(children, parents)]
	for j in desc
		if j > N
			STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab, ind)
		else
			k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
	end
	k, STRING = cp(")", k, STRING)
	k, STRING = cp(nodelab[1], k, STRING)
	k, STRING = cp(";", k, STRING)
	outstring = ""
	for i in STRING
		outstring = "$outstring$i"
	end
	if name == true && prefix != ""
		namebit = chop(match(r"^[^\(]+\(", outstring).match)
		replace(outstring, namebit, "$namebit ")
	end
	return outstring
end
