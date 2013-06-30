# Function for reading in a tree from file.
function TreeRead(filepath::ASCIIString)
	instream = open(expanduser(filepath))
	instring = readall(instream)
	close(instream)
	if search(instring, ";") != 0:-1 && search(instring, "(") != 0:-1 && search(instring, "(")
		trees = split(instring, ';')
		trees = [replace(i, r"(\r|\n|\s)", "") for i in trees]
		trees = trees[bool([t != "" for t in trees])]
		if length(trees) == 1
			if search(trees[1], ":") == 0:-1
				tree = CladoBuild(trees[1])
				return tree
			elseif search(trees[1], ":") != 0:-1
				tree = TreeBuild(trees[1])
				return tree
			end
		end
		outputTrees = Array(Phylogeny, length(trees))
		for i in 1:length(trees)
			if search(trees[i], ":") == 0:-1
				outputTrees[i] = CladoBuild(trees[i])
			elseif search(trees[i], ":") != 0:-1
				outputTrees[i] = TreeBuild(trees[i])
			end
		end
		return outputTrees
	elseif search(instring, "phyloxml") != 0:-1 && search(instring,"<") != 0:-1 && search(instring, ">")
		instring = replace(instring, r"<phyloxml.*>", "")
		trees = split(instring, "</phylogeny>")
		trees = [replace(i, r"(\r|\n|\s)", "") for i in trees]
		trees = trees[bool([t != "" for t in trees])]
		if length(trees) == 1
			tree = phyxbuild(trees[1])
			return tree
		end
		outputTrees = Array(PhyloX, length(trees))
		for i in 1:length(trees)
			outputTrees[i] = phyxbuild(trees[i])
		end
		return outputTrees
	end
end


function phyxbuild(input::ASCIIString)
	inputArray = split(input, r"><")
	inputArray = ["<$i>" for i in inputArray]
end






# Sub function for creation of a Clado structs from newick format strings.
# Used to build Clado structs from newick strings during operation of the 
# TreeRead function.
function CladoBuild(tp::ASCIIString)
	function AddInternal(edge, currentNode, node, index, j)
		edge[j, 1] = currentNode
		node += 1
		edge[j, 2] = currentNode = node
		index[node] = j
		j += 1
		return currentNode, node, j
	end
	function AddTerminal(edge, currentNode, tip, index, tipLabel, tpc, k, j)
		edge[j, 1] = currentNode
		edge[j, 2] = tip
		index[tip] = j
		tipLabel[tip] = tpc[k]
		k += 1
		tip += 1
		j += 1
		return currentNode, tip, k, j
	end
	function GoDown(index, currentNode, nodeLabel, nbTip, tpc, k, edge)
		l = index[currentNode]
		nodeLabel[currentNode - nbTip] = tpc[k]
		k += 1
		currentNode = edge[l, 1]
		return k, currentNode
	end
	tp = "$tp;"
	if ismatch(r"^[^\(]+\(", tp)
		cutoff = length(match(r"^[^\(]+\(", tp).match)
		treeName = chop(match(r"^[^\(]+\(", tp).match)
		treeName = treeName[length(treeName)] == ' ' ? treeName[1:length(treeName)-1] : treeName
		tp = tp[cutoff:length(tp)]
	else treeName = ""
	end
	if search(tp, ",") == 0:-1
		edge = Array(Int, 2,2)
		edge[1,1:2] = [2,1]
		edge[2,1:2] = [1,2]
		tp = split(tp, r"[\\(\\):;]")
		edgeLength = tp[3]
		Nnode = 1
		tipLabel = tp[2]
		if tp[4] != ""
			nodeLabel = tp[4]
		end
		phyloobject = Phylo(edge, Nnode, tipLabel, edgeLength, nodeLabel, -1.0)
		return phyloobject
	end
	tsp = split(tp, "")
	tp = replace(tp, "\s", "")
	tp = replace(tp, ")", ")NA")
	tp = replace(tp, "(", "rem(")		
	tpc = split(tp, r"[\\(\\),;]")
	tpc = tpc[1:length(tpc)-1]
	tpc = tpc[tpc .!= "rem"]
	skeleton = tsp[bool([i == "(" || i == ")" || i == "," || i == ";" for i in tsp])]
	nsk = length(skeleton) # Length of the nexus skeleton.
	nbNode = sum(skeleton .== ")")
	nbTip = sum(skeleton .== ",") + 1
	nbEdge = nbNode + nbTip # Number of edges for tree is number of tips and nodes.
	nodeLabel = ["" for i in 1:nbNode]
	tipLabel = ["" for i in 1:nbTip]
	edge = Array(Int, nbEdge,2)
	currentNode = node = nbTip + 1
	edge[nbEdge, 1] = 0
	edge[nbEdge, 2] = node
	index = [0 for i in 1:nbEdge+1]
	index[node] = nbEdge
	j = k = tip = 1
	for i in 2:nsk
		if skeleton[i] == "("
			currentNode, node, j = AddInternal(edge, currentNode, node, index, j)
		end
		if skeleton[i] == "," && skeleton[i-1] != ")"
			currentNode, tip, k, j = AddTerminal(edge, currentNode, tip, index, tipLabel, tpc, k, j)
		end
		if skeleton[i] == ")"
			if skeleton[i - 1] == ","
				currentNode, tip, k, j = AddTerminal(edge, currentNode, tip, index, tipLabel, tpc, k, j)
				k, currentNode = GoDown(index, currentNode, nodeLabel, nbTip, tpc, k, edge)
			end
			if skeleton[i - 1] == ")"
				k, currentNode = GoDown(index, currentNode, nodeLabel, nbTip, tpc, k, edge)
			end
		end
	end
	edge = edge[1:nbEdge-1, 1:2]
	nodeLabel = [replace(i, r"^NA", "") for i in nodeLabel]
	phyloobject = Clado(treeName, edge, tipLabel, nbNode, nodeLabel)
	return phyloobject
end



# Sub function for creation of a Phylo structs from newick format strings.
# Used to build Phylo structs from newick strings during operation of the 
# TreeRead function.
function TreeBuild(tp::ASCIIString)
	function AddInternal(edge, j, currentNode, node, index)
		edge[j, 1] = currentNode
		node += 1
        edge[j, 2] = currentNode = node
        index[node] = j
        j += 1
        return node, currentNode, j
	end
	function AddTerminal(edge, j, currentNode, tip, index, tpc, k, tipLabel, edgeLength)
		edge[j, 1] = currentNode
        edge[j, 2] = tip
        index[tip] = j
        X = split(tpc[k], ":")
        tipLabel[tip] = X[1]
        edgeLength[j] = X[2]
        k += 1
        tip +=1
        j += 1
        return k, tip, j
	end
	function GoDown(index, currentNode, tpc, k, nodeLabel, nbTip, edgeLength)
		l = index[currentNode]
        X = split(tpc[k], ":")
        nodeLabel[currentNode - nbTip] = X[1]
        edgeLength[l] = X[2]
        k += 1
        currentNode = edge[l, 1]
        return currentNode, k
	end
	tp = "$tp;"
	if ismatch(r"^[^\(]+\(", tp)
		cutoff = length(match(r"^[^\(]+\(", tp).match)
		treeName = chop(match(r"^[^\(]+\(", tp).match)
		treeName = treeName[length(treeName)] == ' ' ? treeName[1:length(treeName)-1] : treeName
		tp = tp[cutoff:length(tp)]
	else treeName = ""
	end
	if search(tp, ",") == 0:-1
		edge = Array(Int, 2,2)
		edge[1,1:2] = [2,1]
		edge[2,1:2] = [1,2]
		tp = split(tp, r"[\\(\\):;]")
		edgeLength = tp[3]
		Nnode = 1
		tipLabel = tp[2]
		if tp[4] != ""
			nodeLabel = tp[4]
		end
		phyloobject = Phylo(edge, Nnode, tipLabel, edgeLength, nodeLabel, -1.0)
		return phyloobject
	end
	tsp = split(tp, "")
	if ismatch(r"\)[^\(:\)\r\n]*;", tp) # Match the root if there is no colon.
		m = match(r"\)[^\(:\)\r\n]*;", tp)
		mstring = m.match[1:length(m.match)-1]
		newString = "$mstring:NA;"
		st1 = tp[1:m.offset-1]
		tp = "$st1$newString"
	end
	tp = replace(tp, ")", ")NA")
	tp = replace(tp, "\s", "")
	tp = replace(tp, "(", "rem(")
	tpc = split(tp, r"[\\(\\),;]")
	tpc = tpc[1:length(tpc)-1]
	tpc = tpc[tpc .!= "rem"]
    skeleton = tsp[bool([i == "(" || i == ")" || i == "," || i == ";" for i in tsp])]
    nsk = length(skeleton)
    nbNode = sum(skeleton .== ")")
    nbTip = sum(skeleton .== ",") + 1
    nbEdge = nbNode + nbTip
    nodeLabel = ["" for i in 1:nbNode]
	tipLabel = ["" for i in 1:nbTip]
	edgeLength = ["" for i in 1:nbEdge]
    edge = Array(Int, nbEdge, 2)
    currentNode = node = nbTip + 1
    edge[nbEdge, 1] = 0
    edge[nbEdge, 2] = node
    index = [0 for i in 1:(nbEdge + 1)]
    index[node] = nbEdge
    j = k = tip = 1
	for i in 2:nsk
        if skeleton[i] == "(" 
            node, currentNode, j = AddInternal(edge, j, currentNode, node, index)
        end
        if skeleton[i] == ","
            if skeleton[i - 1] != ")" 
                k, tip, j = AddTerminal(edge, j, currentNode, tip, index, tpc, k, tipLabel, edgeLength)
            end
        end
        if skeleton[i] == ")"
            if skeleton[i - 1] == ","
                k, tip, j = AddTerminal(edge, j, currentNode, tip, index, tpc, k, tipLabel, edgeLength)
                currentNode, k = GoDown(index, currentNode, tpc, k, nodeLabel, nbTip, edgeLength)
            end
            if skeleton[i - 1] == ")"
                currentNode, k = GoDown(index, currentNode, tpc, k, nodeLabel, nbTip, edgeLength)
            end
        end
    end
	edge = edge[1:nbEdge-1, 1:2]
	rootEdge = edgeLength[nbEdge]
	if rootEdge != "NA"	# Resolve whether there is a rootedge for the tree.
		rootEdge = float64(rootEdge)
	else
		rootEdge = -1.0
	end
	edgeLength = float64([i == "" ? -1.0 : float64(i) for i in edgeLength[1:nbEdge-1]])
	nodeLabel = [replace(i, r"^NA", "") for i in nodeLabel]
	phyloobject = Phylo(treeName, edge, nbNode, tipLabel, edgeLength, nodeLabel, rootEdge)
	return phyloobject
end



# Function method for writing a single Phylogeny type to file.
function TreeWrite(tree::Phylogeny, file::ASCIIString = "output.nwk", append::Bool = false, treeNames::Bool = true)
	output = writenewick(tree, treeNames)
	if append == true
		outstream = open(file, "a")
	else
		outstream = open(file, "w")
	end
	println(outstream, output)
	close(outstream)
end


# Function method for writing multiple trees to file. In order to do this they must be as an array.
function TreeWrite(tree::Array{Phylogeny}, file::ASCIIString = "output.nwk", append::Bool = false, treeNames::Bool = false)
	outputarray = Array(ASCIIString, length(tree))
	for i in 1:length(tree)
		outputarray[i] = writenewick(tree[i], treeNames)
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



# function used by the function which write newick trees to check labels.
function checklabels(labels)
	labels = replace(labels, r"^[[:space:]\\(]+", "")
	labels = replace(labels, r"[[:space:]\\)]+$", "")
	labels = replace(labels, r"[[:space:]]", "_")
	labels = replace(labels, r"[,:;]", "")
	labels = replace(labels, r"[\\(\\)]", "-")
	labels = replace(labels, r"_{2,}", "_")
	labels = replace(labels, r"-{2,}", "-")
	return labels
end



# Function used by the functions which write newick trees.
function cp(x, k, STRING)
	STRING[k] = x
	k += 1
	return k, STRING
end



# Returns the internal node which is the root of the tree.
function getRoot(children, parents)
	r = parents[[find(contains(children, i)) == [1] ? false : true for i in parents]][1]
	return r
end



# This function returns an array of arrays, showing which children
# nodes have.
function getKids(phy::Phylogeny)
	N = length(phy.tipLabel)
	kids = Array(Array{Int64}, N + phy.Nnode)
	for i in 1:N + phy.Nnode
		logic = bool([n == i for n in phy.edge[1:size(phy.edge,1), 1]])
		kids[i] = phy.edge[logic,2]
	end
	return kids
end



function writenewick(phy::Clado, name)
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
	kids = getKids(phy)
	LS = (4 * N) + 5
	LS = LS + N # if there are nodelabels.
	ind = [findin(children, i) for i in 1:max(phy.edge)]
	STRING = ["" for i in 1:LS]
	k = 1
	k, STRING = cp(prefix, k, STRING)
	k, STRING = cp("(", k, STRING)
	root = getRoot(children, parents)
	desc = kids[root]
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
function writenewick(phy::Phylo, name)
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
	kids = getKids(phy)
	LS = (4 * N) + 5
	LS = LS + N
	LS = LS +(4 * N)
	ind = [findin(children, i) for i in 1:max(phy.edge)]
	STRING = ["" for i in 1:LS]
	k = 1
	k, STRING = cp(prefix, k, STRING)
	k, STRING = cp("(", k, STRING)
	root = getRoot(children, parents)
	desc = kids[root]
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