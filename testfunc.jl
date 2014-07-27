function testfunc()
a::PhyNode = PhyNode("A.")
b::PhyNode = PhyNode("B.")
c::PhyNode = PhyNode("C.")
d::PhyNode = PhyNode("D.")
e::PhyNode = PhyNode("E.")
graft!(e, d)
graft!(e, c)
graft!(d, a)
graft!(d, b)
myTree = Phylogeny("My Simple Tree", e, true, false)
println("Node a is a leaf? $(isLeaf(a))")
println("Node b is a leaf? $(isLeaf(b))")
println("Node c is a leaf? $(isLeaf(c))")
println("Node d is a leaf? $(isLeaf(d))")
println("Node e is a leaf? $(isLeaf(e))")
println("Node a is a root? $(isRoot(a))")
println("Node b is a root? $(isRoot(b))")
println("Node c is a root? $(isRoot(c))")
println("Node d is a root? $(isRoot(d))")
println("Node e is a root? $(isRoot(e))")
println("Node a is a node? $(isNode(a))")
println("Node b is a node? $(isNode(b))")
println("Node c is a node? $(isNode(c))")
println("Node d is a node? $(isNode(d))")
println("Node e is a node? $(isNode(e))")
end


                     |======= b
              =======d
              |      |======= a
e =============
              |
              =============== c


function testfunc2()
  a = PhyNode("A.")
  b = PhyNode("B.")
  c = PhyNode("C.")
  d = PhyNode("D.")
  e = PhyNode("E.")
  graft!(e, d)
  graft!(e, c)
  graft!(d, a)
  graft!(d, b)
  myTree = Phylogeny("My Simple Tree", e, true, false)
  traverser = BreadthFirstTraverser(myTree)
  println("Looking at node $(getName(getCurrent(traverser)))")
  while !hasReachedEnd(traverser)
    next!(traverser)
    println("Looking at node $(getName(getCurrent(traverser)))")
  end
  searchDF(myTree, isLeaf)
end