# Dependency for XML.
using LightXML


## Some PhyX 'Extensions' - In this example, I want to be able to read in the tree structure from an XML file, and taxonomy info.
## So I begin by defining a Julia type to store this info and a function that will construct this from xml.
type ID
  provider::String
  identifier::String
end

type Taxonomy
  IDs::Array{ID}
  Code::String
  ScientificName::String
end

function readPhyXML(xml::XMLElement, ::Type{Taxonomy})
  taxxml = get_elements_by_tagname(xml, "taxonomy")
  if !isempty(taxxml)
    idxml = get_elements_by_tagname(taxxml[1], "id")
    if !isempty(idxml)
      idarray = [ID(attribute(i, "provider"; required=false), content(i)) for i in idxml]
    else
      idarray = Array(ID, 0)
    end
    code = get_elements_by_tagname(taxxml[1], "code")
    if !isempty(code)
      codeval = content(code[1])
    else
      codeval = ""
    end
    sciname = get_elements_by_tagname(taxxml[1], "scientific_name")
    if !isempty(sciname)
      name = content(sciname[1])
    else
      name = ""
    end
    return Taxonomy(idarray, codeval, name)
  else
    return Taxonomy(Array(ID, 0), "", "")
  end
end


## This is the core PhyX stuff, the stuff that will get extended with the code above, that extends it's capability to read taxonomy info from Phyloxml files.

# Macro is used to extend the definition of PhyXElement.
macro ExtElement(attributes...)
  code = :(begin end)
  for i in 1:length(attributes)
    push!(code.args, :($(symbol("$(attributes[i])"))::$(attributes[i])))
  end
  code
end

# Nodes of a tree, at most basic each has a label, wether it is the root, a tip, and which node is its parent.
type PhyXElement
  Label::String
  Root::Bool
  Tip::Bool
  Parent::Int
  @ExtElement Taxonomy

  PhyXElement(lab, isroot, istip, parent) = new(lab, isroot, istip, parent) # Constructor allows for incomplete init.
end

# Simple type for storing a tree.
type PhyXTree
  name::String
  nodes::Array{PhyXElement}
  isRooted::Bool
end

# Type used in tree building functions.
type nodeTracker
  nodeIndex::Int
end

macro gatherPhyXMLAttributes()
  code = :(begin end)
  for i in names(PhyXElement)[5:length(names(PhyXElement))]
    push!(code.args, :(cladeArray[currentClade.nodeIndex].$i = readPhyXML(xmlclade, $(i))))
  end
  code
end

function recursiveBuild(xmlclade, cladeArray, currentClade, parentClade::Int)
  # Update the node tracker.
  currentClade.nodeIndex += 1
  current = currentClade.nodeIndex # Initialize a local variable called current, taken from the currentClade variable to keep as the variable to pass to furthur recursive calls as the parent index.
  # Get name of clade element.
  label = ""
  isroot = parentClade == 0 ? true : false
  children = get_elements_by_tagname(xmlclade, "clade")
  istip = length(children) == 0 ? true : false
  # Incompletely initialize the clade object:
  cladeArray[currentClade.nodeIndex] = PhyXElement(label, isroot, istip, parentClade)
  # Get and process all additional data and complete the clade object.
  @gatherPhyXMLAttributes
  # Make the clade object by correctly calling the constructor.
  for i in children
    recursiveBuild(i, cladeArray, currentClade, current)
  end
end

function readPhyXML(xmltree::XMLElement, ::Type{PhyXTree})
  treestring = replace(string(xmltree), r"(\r|\n|\t)", "")
  treestring = replace(treestring, r"(\s{2,})", "")
  tstable = split(treestring, "><")
  startclade = 0
  endclade = 0
  for i in tstable
    if i == "clade"
      startclade += 1
    end
    if i == "/clade"
      endclade += 1
    end
  end
  startclade != endclade ? println("Warning! There are unequal numbers of clade begins and clade ends in this tree") : println("Current tree has $startclade clade nodes")
  # Ok the number of nodes has been established.
  Clade = Array(PhyXElement, startclade)    # Make an array to contain the Clade elements.
  BackTrack = zeros(Int, startclade)      # Make an array which tracks the parent of a Clade, to allow backtracking.
  Current = nodeTracker(0)                    # Start the nodetracker type.
  BackTrack[1] = 0
  XML = get_elements_by_tagname(xmltree, "clade")[1]
  extensions = names(PhyXElement)[5:length(names(PhyXElement))]
  recursiveBuild(XML, Clade, Current, 0)
  treename = ""
  isrooted = false
  return PhyXTree(treename, Clade, isrooted)
end

function readtree(file::ASCIIString, ::Type{PhyXTree})
  treedoc = parse_file(expanduser(file))
  phylogenies = get_elements_by_tagname(root(treedoc), "phylogeny")
  if length(phylogenies) > 1
    return PhyXTree[readPhyXML(i, PhyXTree) for i in phylogenies]
  else
    return readPhyXML(phylogenies[1], PhyXTree)
  end
  println("No Trees found in file...")
end



