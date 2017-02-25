xxxxa=whos;
%xxxxnames=a.name;
%xxxxclasses=a.class;

for xxxxk=1:length(xxxxa)
    if strcmp(xxxxa(xxxxk).class,'double')
        
        assignin('base','xxxtest',b(:,idx2(k)));
        
        xxxxa(xxxxk).name
        s=eval('unique(xxxxa(xxxxk).name)')
        pause
        if length(s)==2            
        if s(2)==1&&s(1)==0
            eval('xxxxa(xxxxk).name=logical(xxxxa(xxxxk).name)')            
        end            
        end
    end
end
clear xxxa xxxxk;