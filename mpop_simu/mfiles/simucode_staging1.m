function [segseq,segsites]=simucode
N=200;  L=1000;
u=1e-2; r=1e-2;
s=0.01; h=0.5;
w=[1 1+s*h 1+s];
w=w./max(w);

HGa=logical(sparse(2*N,L));
HGb=logical(sparse(2*N,L));

for G=1:200
    
    % mutation
    mutnum=poissrnd(u*N);
    for k=1:mutnum
        i=ceil(rand*2*N);
        j=ceil(rand*L);
        HGa(i,j)=true;
    end
    
    c=0;
    %while c<N
    for c=1:4:2*N
	    s1=ceil(rand*N); s2=ceil(rand*N);
	    while(s2==s1), s2=ceil(rand*N); end
	    g1=s1*2-1; g2=s1*2; g3=s2*2-1; g4=s2*2;
        
        % recombination [hap1,hap2]=HGa(s1,s2);
        gm1=HGa(g1,:); gm2=HGa(g2,:);
        gm3=HGa(g3,:); gm4=HGa(g4,:);
        recombnum1=poissrnd(r*N);
        recombnum2=poissrnd(r*N);
        for k=1:recombnum1
            j=ceil(rand*L);
            gm1new=[gm1(1,1:j),gm2(1,j+1:end)];
            gm2new=[gm2(1,1:j),gm1(1,j+1:end)];
            gm1=gm1new;
            gm2=gm2new;
        end
        for k=1:recombnum2
            j=ceil(rand*L);
            gm3new=[gm3(1,1:j),gm4(1,j+1:end)];
            gm4new=[gm4(1,1:j),gm3(1,j+1:end)];
            gm3=gm3new;
            gm4=gm4new;
        end
	    if (rand>0.5), child1gm1=gm1; else child1gm1=gm2; end
	    if (rand>0.5), child1gm2=gm3; else child1gm2=gm4; end
	    if (rand>0.5), child2gm1=gm1; else child2gm1=gm2; end
	    if (rand>0.5), child2gm2=gm3; else child2gm2=gm4; end

        w(1+child1gm1(selsite)+child1gm2(selsite))
        w(1+child2gm1(selsite)+child2gm2(selsite))

	    if (rand>0.5), HGb(c,:)=gm1; else HGb(c,:)=gm2; end
	    if (rand>0.5), HGb(c+1,:)=gm3; else HGb(c+1,:)=gm4; end
	    if (rand>0.5), HGb(c+2,:)=gm1; else HGb(c+2,:)=gm2; end
	    if (rand>0.5), HGb(c+3,:)=gm3; else HGb(c+3,:)=gm4; end
    end
    
    HGa=HGb;
    spy(HGa,'*');
    F(G) = getframe;
end

x=sum(HGa);
segsites=find(x>0);
segseq=full(HGa(:,segsites));

%%
% movie(F)
%{
x=sum(HGa);
segsites=find(x>0);
segseq=full(HGa(:,segsites));
freqv=full(x(segsites))./(2*N);

figure;
histsfs(freqv,2);

fprintf('//\nsegsites: ');
fprintf('%d ',segsites);
fprintf('\n');
for k=1:size(segseq,1)
    fprintf('%d',full(segseq(k,:)));
    fprintf('\n');
end
%}

%%

%simusum(ms_mex(400,0,19),N)
%tajima89d_test(ms_mex(400,0,19)+1)
%d=tajima89d_test(ms_mex(400,0,19)+1)


