## Types for tree annotations...
type ID
  Provider::String
  Value::String
end

type URI
  Description::String
  Type::String
  Identifier::String
end

type Taxonomy
  IdSource
  Id::ID
  Code
  ScientificName::String
  Authority::String
  CommonName::String
  Synonym::String
  Rank::String
  Uri::URI
end

## Recurvive extendible type for representation of phylogenetic trees in in Julia.

# Extension type - parametric, contains
type PhyXExtension{T}
  value::T
end

type PhyXElement
  Label::String
  Root::Bool
  BranchLength::Float64
  Extensions::Array{PhyXExtension, 1}
  Children::Array{PhyXElement, 1}
  Parent::PhyXElement
  PhyXElement(label::String, root::Bool, branchlength::Float64, ext::Array{PhyXExtension, 1}) = new(label, root, branchlength, ext, PhyXElement[])
end

function PhyXElement(label::String, root::Bool, branchlength::Float64, ext::Array{PhyXExtension, 1}, parent::PhyXElement)
  x = PhyXElement(label, root, branchlength, ext)
  x.Parent = parent
  return x
end


## Node Manipulation / methods on the PhyXElement type...

# Label manipulation...
function getLabel(x::PhyXElement)
  return x.Label
end

function setLabel!(x::PhyXElement, label::String)
  x.Label = label
end

function getParent(x::PhyXElement)
  return x.Parent
end

# Find out if  node has children, return an array of refs to the children of the node, add children, and remove children of a given node.
function hasChildren(x::PhyXElement)
  return length(x.Children) > 0 ? true : false
end

function getChildren(x::PhyXElement)
  return x.Children
end

function addChild!(parent::PhyXElement, child::PhyXElement)
  if !in(child, parent.Children)
    push!(parent.Children, child)
    if isdefined(child, :Parent)
      child.Parent.Children = child.Parent.Children[child.Parent.Children .!= child]
    end
    child.Parent = parent
  end
end

function removeChild!(Parent::PhyXElement, Child::PhyXElement)
  filter!(x -> x == Child, Parent.Children)
end

function removeChild!(Parent::PhyXElement)
  Parent.Children = PhyXElement[]
end

# Find if given node is a Leaf node, or is the root node, set a node as the root.
function isLeaf(x::PhyXElement)
  return !hasChildren(x)
end

function isRoot(x::PhyXElement)
  return x.Root
end

# Needs to be improved so as the tree is rearranged correctly.
function setAsRoot!(x::PhyXElement)
  x.Root = true
end

function showSiblings(x::PhyXElement)
  return getChildren(x.Parent)
end

function addSiblings!(x::PhyXElement, Siblings::PhyXElement...)
    for i in Siblings
      addChild!(x.Parent, i)
    end
end


## Tree type.
type PhyXTree
  Name::String
  Root::PhyXElement
  Rooted::Bool
  Rerootable::Bool
end


## Methods on the Tree type:

## Depth first tree search.
# Depth first search function uses pre-order depth first search. If no node satisfies the conditions then Nothing is returned.
# Search procedures are iterative for both DF and BF searches, largely the only difference between them is DF uses a Stack (FILO) and
# breadth first uses a Queue (FIFO).
# The conditions variable is a tuple of functions.
# For each type of search BF and DF, there is a searchAny and a searchAll function. The searchAny functions return references to the nodes that returned true for at least one
# of the functions provided in Conditions..., searchAll


function searchDF(Phylogeny::PhyXTree, Condition::Function)
  stack::Stack{PhyXElement} = Stack(PhyXElement)
  push!(stack, Phylogeny.Root)
  while length(stack) > 0
    current::PhyXElement = pop!(stack)
    if Condition(current)
      return current
    else
      for i in current.Children
        push!(stack, i)
      end
    end
  end
end

function searchAllDF(Phylogeny::PhyXTree, Condition::Function)
  stack::Stack = Stack(PhyXElement)
  matches::Array{PhyXElement, 1} = PhyXElement[]
  push!(stack, Phylogeny.Root)
  while length(stack) > 0
    current::PhyXElement = pop!(stack)
    if Condition(current)
      push!(matches, current)
    end
    for i in current.Children
      push!(stack, current)
    end
  end
  return matches
end

function searchBF(Phylogeny::PhyXTree, Condition::Function)
  queue::Queue{PhyXElement} = Queue(PhyXElement)
  enqueue!(queue, Phylogeny.Root)
  while length(queue) != 0
    current::PhyXElement = dequeue!(queue)
    if any([i(current) for i in Conditions])
      return current
    else
      for i in current.Children
        enqueue!(queue, i)
      end
    end
  end
end

function searchAllBF(Phylogeny::PhyXTree, Condition::Function)
  queue::Queue{PhyXElement} = Queue(PhyXElement)
  matches::Array{PhyXElement, 1} = PhyXElement[]
  enqueue!(queue, Phylogeny.Root)
  while length(queue) != 0
    current::PhyXElement = dequeue!(queue)
    if all([i(current) for i in Conditions])
      push!(matches, current)
    else
      for i in current.Children
        enqueue!(queue, i)
      end
    end
  end
  return matches
end
