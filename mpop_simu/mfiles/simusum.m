function simusum(H1,N)

if issparse(H1)
    x=sum(H1);
    segsites=find(x>0);
    segseq=full(H1(:,segsites));
    freqv=full(x(segsites))./(2*N);
else
    segseq=logical(H1);
    freqv=sum(H1)./(2*N);
end

figure;
histsfs(freqv,2);

fprintf('//\nsegsites: ');
%fprintf('%d ',segsites);
fprintf('\n');
for k=1:size(segseq,1)
    fprintf('%d',full(segseq(k,:)));
    fprintf('\n');
end