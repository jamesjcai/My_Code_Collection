function [idx] = inladder2(ladder,posx,moreone)
%[idx] = inladder2(ladder,posx,moreone)

if nargin<3
    moreone=true;
end

idx=[];

for k=1:length(posx)

    pos=posx(k);
id=find(pos==ladder(:,1));
if isempty(id)
id=find(pos==ladder(:,2));
if isempty(id)    
    x=(pos>=ladder);
    id=find(x(:,1)*10+x(:,2)==10);
     if(isempty(id))
         id=0;
     end
end
end
if moreone
    idx=[idx;id];
else
    idx=[idx,id(1)];
end
end
