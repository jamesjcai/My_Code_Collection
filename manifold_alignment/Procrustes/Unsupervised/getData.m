function [X1, X2, X3]=getData
    file1='1G7O-1.pdb.3d.real';
    X1=readProtein(file1);
    meanX1=mean(X1);
    X1=X1-repmat(meanX1, 215,1);
    
    file2='1G7O-21.pdb.3d.real';
    X2=readProtein(file2); X2=4*X2;
    meanX2=mean(X2);
    X2=X2-repmat(meanX2, 215,1);
    
    file3='1G7O-10.pdb.3d.real';
    X3=readProtein(file3); X3=2*X3;
    meanX3=mean(X3);
    X3=X3-repmat(meanX3, 215,1);
    X1=X1';
    X2=X2';
    X3=X3';
end

function X=readProtein(file_name)
    %read a protein
    Fid=fopen (file_name, 'r');
    data=fscanf(Fid, '%f');
    dim = data(1);
    B=data(2:dim*3+1);
    X=(reshape(B,3,dim))';
    fclose(Fid);
end