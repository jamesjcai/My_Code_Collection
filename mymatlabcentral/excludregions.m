function res=excludregions(regv,startn,endn)

%region vector including startn and endn
%for example
%regv=[22 30; 40 50];
%maxn=100;
%res=excludregions(regv,1,100)
%res=    1    21
%    31    39
%    51   100

regv=unique(regv,'rows');

startn=max(startn,1);
res=[startn, regv(1,1)-1];
for (k=1:size(regv,1)-1)
	a=regv(k,2)+1;
	b=regv(k+1,1)-1;
	if (b>a)
		res=[res; [a, b]];
	end
end
res=[res; [regv(end,2)+1, endn]];
