function [basneib]=snpancestor(matfileid,pos)

len=length(pos);
%bas=5*ones(1,len,1);    % ACGTN, 6=N
basneib=6*ones(len,3);    % ACGTN, 6=N


%matfile=['y:/chimpaln/1Mbmat/aln',int2str(matfileid),'.mat'];
matfile=['C:/biodata/chimpaln/1Mbmat_chimp/aln',int2str(matfileid),'.mat'];
if ~(exist(matfile,'file')),
    return    
%error('nofile')
end

eval(['load ',matfile,' xpos xseq'])

%system(['copy ',matfile,' .\\xxx.mat']);
%load xxx xpos xseq
%delete('xxx.mat')

for (k=1:len),
   p=pos(k);                         % this SNP position such as 123456
   idx=find(sum(xpos>p,2)==1);
   errflag=0;

   if ~(isempty(idx)),

       seq=xseq{idx};
       seq(:,find(seq(1,:)==nt2int('-')))=[];   % remove human gaps
       loc=p-xpos(idx,1)+1;           % relative location of site   123456 - 120000 + 1
       
       a=xpos(idx,2)-xpos(idx,1)+1;
       b=size(seq,2);
       if (abs(a-b)>0),
          fprintf(['%d~=%d, Something wrong at %d position %d\n'],a,b,matfileid,p);
	  %error('here')
	  errflag=1;    
       end

       if (loc>1)
	       basneib(k,1)=seq(2,loc-1);
	       basneib(k,2)=seq(2,loc);
	       basneib(k,3)=seq(2,loc+1);
       end
   end
end

%disp(sprintf(['Something wrong at %d\n'],matfileid));

