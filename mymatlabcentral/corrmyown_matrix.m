function [R,P]=corrmyown_matrix(listvar)
%[R,P]=corrmyown_matrix(listvar)
if ~iscell(listvar)
    a=listvar;
    b{size(a,2)}=a(:,end);
    for k=1:size(a,2)-1
        b{k}=a(:,k);
    end
    listvar=b;
end



n=length(listvar);
R=zeros(n);
P=zeros(n);

for i=1:n-1
for j=i+1:n
    [a,b]=corrmyown(listvar{i},listvar{j},[],true);
    R(i,j)=a;
    P(i,j)=b;
    R(j,i)=a;
    P(j,i)=b;    
end
end

%xx={loso,con,cln,cbc,cif,cnm,cco,dN,dS,dN./dS};
%rlab='loso con cln cbc cif cnm cco dN dS dNdS';
%printmat(R,'name',rlab,rlab)
%name={'loso','con','cln','cbc','cif','cnm','cco','dN','dS','dNdS'};


for i=1:n
for j=1:n    
    fprintf('%f\t',R(i,j));
end
    fprintf('\n');
end

fprintf('\n');

for i=1:n
for j=1:n    
    fprintf('%d\t',P(i,j));
end
    fprintf('\n');
end