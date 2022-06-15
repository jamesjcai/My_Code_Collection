delete("my_example_filex.h5")
%delete("my_example_file.h5")

%%
testdat1 = rand(10)*10;
h5create('my_example_filex.h5', '/c', size(testdat1));
h5write('my_example_filex.h5', '/c', testdat1);

testdat2=["geneA";"geneB"];
h5create('my_example_filex.h5', '/s', size(testdat2),'Datatype','string');
h5write('my_example_filex.h5', '/s', testdat2);

save aaa testdat1 testdat2 -v7.3

writematrix(testdat2,'genelist.txt');



