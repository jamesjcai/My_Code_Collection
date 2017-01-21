function [segseq,segsites,freqx]=simucode_noL
N=200;  L=1000;
u=1e-2; r=1e-2;
s=0.4; h=0.5;
w=[1 1+s*h 1+s];
w=w./max(w);

%HGa=logical(sparse(2*N,L));
%HGb=logical(sparse(2*N,L));
HGa=false(2*N,L);
HGb=false(2*N,L);
selsite=500;
[a,b]=ms_mex(400,4,0,4,2);
b=ceil(b*L);
HGa(:,b)=a;


totalG=200;
freqx=[];

for G=1:totalG
    % add selected site
    if ~any(HGa(:,selsite))
        i=ceil(rand*2*N);
        HGa(i,selsite)=true;
    end
    if all(HGa(:,selsite))         % fixed
        HGa(:,selsite)=false;
    end
    
    % mutation
    mutnum=poissrnd(u*N);
    for k=1:mutnum
        i=ceil(rand*2*N);
        j=ceil(rand*L);
        HGa(i,j)=true;
    end
    
    c=1;
    while c<=N
	    s1=ceil(rand*N); s2=ceil(rand*N);
	    while(s2==s1), s2=ceil(rand*N); end
	    g1=s1*2-1; g2=s1*2; g3=s2*2-1; g4=s2*2;
        
        % recombination [hap1,hap2]=HGa(s1,s2);
        p1gm1=HGa(g1,:); p1gm2=HGa(g2,:);
        p2gm1=HGa(g3,:); p2gm2=HGa(g4,:);
        recombnum1=poissrnd(r*N);
        recombnum2=poissrnd(r*N);
        for k=1:recombnum1
            j=ceil(rand*L);
            gm1new=[p1gm1(1,1:j),p1gm2(1,j+1:end)];
            gm2new=[p1gm2(1,1:j),p1gm1(1,j+1:end)];
            p1gm1=gm1new;
            p1gm2=gm2new;
        end
        for k=1:recombnum2
            j=ceil(rand*L);
            gm3new=[p2gm1(1,1:j),p2gm2(1,j+1:end)];
            gm4new=[p2gm2(1,1:j),p2gm1(1,j+1:end)];
            p2gm1=gm3new;
            p2gm2=gm4new;
        end
	    if (rand>0.5), childgm1=p1gm1; else childgm1=p1gm2; end
	    if (rand>0.5), childgm2=p2gm1; else childgm2=p2gm2; end

        if w(1+childgm1(selsite)+childgm2(selsite))>rand
            HGb(c*2-1,:)=childgm1;
            HGb(c*2,:)=childgm2;
            c=c+1;
        end
    end
    
    HGa=HGb;
    freqx(G)=sum(HGa(:,selsite))./(2*N);

if nargout<1    
    figure(1)
    subplot(2,1,1)
    spy(HGa,'.');
    %imagesc(HGa);
    title(sprintf('Generation: %d; u=%g; r=%g; s=%g; h=%g',G,u,r,s,h));
    F(G) = getframe;
    
    
    subplot(2,1,2)
    %{
    if sum(HGa(:,selsite))==2*N
        plot(1:G,freqx,'r');
    else
        plot(1:G,freqx,'g');
    end
    xlim([1 totalG]);
    ylim([0 1]);
    %}      
    v=slidingfun(@tajima89d_test,HGa+1,50,10);
    plot(v);
end
end

x=sum(HGa);
segsites=find(x>0);
segseq=full(HGa(:,segsites));

%%
% movie(F)


