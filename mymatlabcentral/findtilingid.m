function id=findtilingid(chrid,startn,endn,X)

%chrid=22;
%startn=48999990;
%endn=49000005;
%chrid=1;
%starn=50232577;
%endn=50234461;

if(nargin<4)
X=tilingschema(500000);
end

idx=X(find(X(:,2)==chrid),1);
coord=X(find(X(:,2)==chrid),[3 4]);

id1=startn>coord;
row1=find(id1(:,1)*10+id1(:,2)==10);
id2=coord>endn;
row2=find(id2(:,1)*10+id2(:,2)==1);

    id=idx([row1,row2]);
    
%if (row1==row2)
%    id=idx(row1,r);
%else
%    id=idx([row1,row2]);
%end
    

