N=200;
L=1000;

u=1e-2;
r=1e-2;

H1=logical(sparse(2*N,L));
H2=logical(sparse(2*N,L));

for G=1:100
    mutnum=poissrnd(u*N);
    for k=1:mutnum
        i=ceil(rand*2*N);
        j=ceil(rand*L);
        H1(i,j)=true;
    end
    
    s1=ceil(rand*N);
    s2=ceil(rand*N);
    while(s2==s1), s2=ceil(rand*N); end
    % [hap1,hap2]=H1(s1,s2);
    g1=s1*2-1; g2=s1*2; g3=s2*2-1; g4=s2*2;
    
    
    % Hap=Hap(randperm(2*N),:);
    % F(G) = getframe;
end
spy(Hap,'*');
% movie(F)
