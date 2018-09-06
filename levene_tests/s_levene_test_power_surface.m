
smpsiz=15:5:115;
lamdav=1:0.1:3;

%smpsiz=5:5:105;
%lamdav=1:0.05:2;

Zc=cell(21,21);
%%
for kx=1:21
    for ky=1:21
        Sz=smpsiz(kx);
        Lv=lamdav(ky);
        p=ones(1000,1);        
        g=[ones(Sz,1);ones(Sz,1)*2];
        parfor k=1:1000
            Ra=randn(Sz,1);
            Rb=randn(Sz,1)*Lv;
            RX=[Ra;Rb];            
            p(k)=vartestn(RX,g,'display','off','testtype','LeveneQuadratic');
        end
        Zc{kx,ky}=p;
    end
end
[X,Y] = meshgrid(smpsiz,lamdav);
% save zxy Zc X Y
%%
Z=cellfun(@(X)(sum(X<0.05)),Zc);
figure;
% mesh(X,Y,Z)
% contour(X,Y,Z./10,'ShowText','on')
[C,h] = contour(X,Y,Z./10,'-k','ShowText','on');
clabel(C,h,'FontSize',12,'Color','k')
%xlabel('Num of samples per group')
%ylabel('Ratio of STD btw groups')
