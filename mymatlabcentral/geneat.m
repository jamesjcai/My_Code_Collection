function [genelist]=geneat(chrid,startn,endn)


%geneat(3,120500001,123400000)
load('y:/humangenes/humangenelist','D');
%[ok,idx]=ismember(id,D.geneid);

[chrid,idx]=sort(chrid);
startn=startn(idx);
endn=endn(idx);

%chrid=[3 4];
%startn=[120500001 120500001];
%endn=[123400000 123400000];


t=1;

uchrid=unique(chrid);

for k=1:length(uchrid)
    
    id=chrid==uchrid(k);    
    chridi=chrid(id);
    startni=startn(id);
    endni=endn(id);
    
   id2=uchrid(k)==D.chrid;
   geneid=D.geneid(id2);
   X=D.transstart(id2);   % c
   Y=D.transend(id2);       % d   
    
    for j=1:length(startni)
        a=startni(j);
        b=endni(j);
        idxx=find((X>=a & X<=b) |  (Y>=a & Y<=b));
        if ~isempty(idxx)
            genelist{t}=geneid(idxx);
        end
        t=t+1;
    end
end
    

t=1;
for k=1:length(genelist)
for j=1:length(genelist{k})
    genelist2{t}=genelist{k}{j};
    t=t+1;
end
end    
genelist=unique(genelist2)';
