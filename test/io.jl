# test @tr_str
begin
    # for Cladograms
    parsedtree = @tr_str  "(A,B,(C,D)E)F;"
    @test parsedtree == Clado("",[5 1
                                  5 2
                                  5 6
                                  6 3
                                  6 4],2,
                              ["A","B","C","D"],
                              ["F","E"])


    # for Phylograms
    parsedtree = @tr_str "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"
    @test parsedtree == Phylo("", [5 1
                                   5 2
                                   5 6
                                   6 3
                                   6 4],2,
                              ["A","B","C","D"],
                              [0.1,0.2,0.5,0.3,0.4],
                              ["F","E"],
                              -1.0)
end

# test writenewick
begin
    # for Phylograms
    @test Phylogenetics.newick(Clado("test",[5 1
                                             5 2
                                             5 6
                                             6 3
                                             6 4],2,
                                     ["A","B","C","D"],
                                     ["F","E"]),false) == "(A,B,(C,D)E)F;"

    # for Cladograms
    @test Phylogenetics.newick(Phylo("test", [5 1
                                              5 2
                                              5 6
                                              6 3
                                              6 4],2,
                                     ["A","B","C","D"],
                                     [0.1,0.2,0.5,0.3,0.4],
                                     ["F","E"],
                                     -1.0),true) == "test(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"
end
