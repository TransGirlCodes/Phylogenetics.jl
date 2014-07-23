# Recursive Build Functions...

# Recusrsiveley builds tree structure from PhyloXML format.
function buildRecursively(xmlclade::XMLElement, extensionsArray, parent::PhyXElement)
  label::String = ""
  children::Array{XMLElement, 1} = get_elements_by_tagname(xmlclade, "clade")
  bl = attribute(xmlclade, "branch_length")
  element::PhyXElement = PhyXElement(label, false, (bl == nothing ? 0.0 : float(bl)), PhyXExtension[parsePhyloXML(xmlclade, i) for i in extensionsArray], parent)
  element.Children = PhyXElement[buildRecursively(n, extensionsArray, element) for n in children] # Have reservations about this line, I reckon to be consistent with the API perhaps addChild!() should be used, however that means iterative uses use push! and I don't know if that's better or worse than this line performance wise.
  return element
end

function buildRecursively(xmlclade::XMLElement, extensionsArray, ::Nothing)
  label::String = ""
  children::Array{XMLElement, 1} = get_elements_by_tagname(xmlclade, "clade")
  element::PhyXElement = PhyXElement(label, true, 0.0, PhyXExtension[parsePhyloXML(xmlclade, i) for i in extensionsArray])
  element.Children = PhyXElement[buildRecursively(n, extensionsArray, element) for n in children]
  return element
end

function parsePhyloXML(xmltree::XMLElement, ::Type{PhyXTree}, extensions) # extensions is a tuple of types.
  # Get the first clade of the tree - the root clade in the XML, the recursiveley build up the tree structure.
  # Function assumes a tree is rerootable unless otherwise specified, and assumes it's not rooted unless otherwise specified in the file.
  # Function also assumes that there is and can only be one name for the tree.
  XML::XMLElement = get_elements_by_tagname(xmltree, "clade")[1]
  root::PhyXElement = buildRecursively(XML, extensions, nothing)
  treename::String = ""
  rerootable::Bool = true
  rooted::Bool = false
  attribute(xmltree, "rooted") == "true" ? rooted = true : nothing
  attribute(xmltree, "rerootable") == "false" ? rerootable = true : nothing
  name::Array{XMLElement, 1} = get_elements_by_tagname(xmltree, "name")
  length(name) == 1 ? treename = content(name[1]) : nothing
  return PhyXTree(treename, root, rooted, rerootable)
end


# Open to suggestions making this function neater re deciding on different formats. Since branching conditionally on the format variable seems less neat to me,
# than taking advantage of method dispatch. Perhaps some dummy types for the input formats to allow method dispatch would be better?

function readPhyloXML(file::String, Extensions...)
  treedoc::XMLDocument = parse_file(expanduser(file))
  phylogenies::Array{XMLElement,1} = get_elements_by_tagname(root(treedoc), "phylogeny")
  if length(phylogenies) > 1
    return [parsePhyloXML(i, PhyXTree, Extensions) for i in phylogenies]
  else
    return parsePhyloXML(phylogenies[1], PhyXTree, Extensions)
  end
  error("No Trees found in file...")
end


# Build up the tree from a newick string by recursive descent.
function buildRecursively(newick::String, from::Int, to::Int, parent::PhyXElement)
  if newick[from] != '\('
    setLabel!(parent, newick[from:to])
  else
    bracketCounter::Int = 0
    colonMarker::Int = 0
    positionMarker::Int = from
    for i in from:to
      token::Char = newick[i]
      if token == '\('
        bracketCounter += 1
      elseif token == '\)'
        bracketCounter -= 1
      elseif token == ':'
        colonMarker = i
      end
      if (bracketCounter == 0) || (bracketCounter == 1 && token == ',')
        addChild!(parent, buildRecursively(newick, positionMarker + 1, colonMarker, PhyXElement("", false, float(newick[(colonMarker + 1):(i - 1)]), PhyXExtension[], parent)))
        positionMarker = i
      end
    end
  end
end



function makeNodeFromNS(SubString::String, Depth::Int)
    if search(SubString, '(') == 0
      error(string("Tree is not well formed, location: ", SubString))
    end
    # Split the Substring into Children.
    childrenStrings = split(SubString, ',')
    childrenNodes = PhyXElement[]
    for cs in childrenStrings
      if length(cs) == 0
        continue
      end
      nodeInformation = split(cs, ':')
      name, bl::Float64, bs::Float64 = "", 0.0, 0.0
      if length(nodeInformation) == 2 # If nodeinfo is of length 2, the node contains info on the name / branchlength, or possibly bootstrap support.
        bl = float(nodeInformation[2])
        name = nodeInformation[1]
        # The nexus format may put bootstrap values in the names position...
        try
          name = float(name)
          if 0 <= name <= 1
            bs = name
          elseif 1 <= name <= 100
            bs = name / 100
          end
          name = ""
        catch
          name = nodeInformation[1]
        end
      else # otherwise assume the length is 1 in which case the entry only contains a name...
        name = nodeInformation[1]
      end
      node = PhyXElement(name, false, bl, PhyXExtension[])
      push!(childrenNodes, node)
    end
    return childrenNodes
end

function parseNewickNode(NewickString, Depth)
  NewickString == "" ? error("Input NewickString is devoid of characters.") : nothing
  search(NewickString, '(') == 0 ? return makeNodesFromNS(NewickString, Depth) : nothing

  nodes = PhyXElement[]
  children = PhyXElement[]
  start::Int = 1
  lengthOfPreNewickString::Int = 0
  bracketStack::Stack = Stack()

  for i in 1:length(NewickString)
    if NewickString[i] == '('
      push!(bracketStack, i)
      continue
    end

    if NewickString[i] == ')'
      j = pop!(bracketStack)
      if length(bracketStack) == 0
        InternalNode = nothing
        subNewick = NewickString[start + lengthOfPreNewickString: j]
        preceedingNodes = makeNodeFromNS(subNewick, Depth)
        nodes = vcat(nodes, preceedingNodes)
        if i + 1 < length(NewickString)
          rightOfBracket = NewickString[i + 1:end]
          m = match(r"[\)\,\(]", rightOfBracket)
          if m != nothing
            indexOfNextSymbol = m.offset
            repOfInternalNode = rightOfBracket[1:indexOfNextSymbol]
            internalNodes = makeNodeFromNS(repOfInternalNode, Depth)
            if length(internalNodes) > 0
              InternalNode = internalNodes[1]
            end
            lengthOfPreNewickString = length(repOfInternalNode)
          else
            InternalNode = makeNodeFromNS(NewickString[i + 1:end], Depth)[1]
            lengthOfPreNewickString = length(NewickString) - i
        end
        if InternalNode == nothing
          InternalNode = PhyXElement("", false, 0.0, PhyXExtension[])
          lengthOfPreNewickString = 0
        end

        childNewickString = NewickString[j + 1 : i]
        addChild!(InternalNode, )
end



function parseNewick(newickString::String, rooted::Bool, rerootable::Bool)
  newickString = replace(newickString, r"(\r|\n|\s)", "")
  treename = ""
  if ismatch(r"^[^\(]+\(", newickString)
    cutoff = length(match(r"^[^\(]+\(", newickString).match)
    treeName = chop(match(r"^[^\(]+\(", newickString).match)
    treeName = treeName[length(treeName)] == ' ' ? treeName[1:length(treeName)-1] : treeName
    newickString = newickString[cutoff:length(newickString)]
  end
  x::Int = rsearchindex(newickString, ":")
  root:: PhyXElement = PhyXElement("", true, 0.0, PhyXExtension[])
  buildRecursively(newickString, 1, x, root)
  return PhyXTree(treename, root, rooted, rerootable)
end






#function buildRecursively(newick::String, parent::PhyXElement)

#end

#function parseNewick(newickString::String, rooted::Bool)

#end


function readNewick(file::String)
  instream = open(expanduser(file))
  instring = readall(instream)
  close(instream)
  newickstrings = split(instring, ';')
  newickstrings = [replace(i, r"(\r|\n|\s)", "") for i in newickstrings]
  newickstrings = newickstrings[bool([length(t) > 0 for t in newickstrings])]
  if length(newickstrings) > 1
    return [parseNewick(i) for i in newickstrings]
  else
    return parseNewick(newickstrings[1])
  end
  error("No phylogenies detected in the file.")
end


















