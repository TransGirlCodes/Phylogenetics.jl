using LibExpat
filepath = "~/Desktop/apaf-1.xml";
instream = open(expanduser(filepath));
instring = readall(instream);
close(instream);
inputstring = instring;
xmltree = xp_parse(inputstring);
allInTrees = xmltree[xpath"phylogeny"];
tree = allInTrees[1];
rooted = get(tree.attr, "rooted", "") == "true" ? true : false;
rerootable = get(tree.attr, "rerootable", "") == "true" ? true : false;
nClades = length(find(tree, "//clade"));
allClades = find(tree, "//clade");
structure = zeros(Int, nClades);