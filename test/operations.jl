# test hashing
begin
    @test hash(Clado("",[1 2],1,["A"],["B"])) ==
          hash(Clado("",[1 2],1,["A"],["B"]))

    @test hash(Phylo("",[1 2],1,["A"],[1],["B"],-1.0)) ==
          hash(Phylo("",[1 2],1,["A"],[1],["B"],-1.0))
end

#test equality
begin
    @test isequal(Clado("",[1 2],1,["A"],["B"]),
                  Clado("",[1 2],1,["A"],["B"]))
    @test isequal(Phylo("",[1 2],1,["A"],[1],["B"],-1.0),
                  Phylo("",[1 2],1,["A"],[1],["B"],-1.0))
end


# test tree querying functions

clado = Clado("test",[5 1
                      5 2
                      5 6
                      6 3
                      6 4],2,
              ["A","B","C","D"],
              ["F","E"])

phylo =Phylo("test", [5 1
                      5 2
                      5 6
                      6 3
                      6 4],2,
             ["A","B","C","D"],
             [0.1,0.2,0.5,0.3,0.4],
             ["F","E"],
             -1.0)


#test getroot
begin
    @test getroot(clado.edge[:,2],clado.edge[:,1]) == 5
    @test getroot(phylo.edge[:,2],phylo.edge[:,1]) == 5
end

#test getkids
begin
    @test getkids(clado) == Array[[], [], [], [],
                                  [1,2,6],
                                  [3,4]]

    @test getkids(phylo) == Array[[], [], [], [],
                                  [1,2,6],
                                  [3,4]]
end
