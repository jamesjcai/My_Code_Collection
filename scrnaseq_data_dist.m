figure;
[m,c]=onekpoiss(1);
hold on
scatter(m,c)
[m,c]=onekpoiss(2.5);
scatter(m,c)
[m,c]=onekpoiss(5);
scatter(m,c)
[m,c]=onekpoiss(10);
scatter(m,c)



figure;
[m,c,h]=onekpoiss(1);
hold on
scatter3(m,c,h)
[m,c,h]=onekpoiss(2.5);
scatter3(m,c,h)
[m,c,h]=onekpoiss(5);
scatter3(m,c,h)
[m,c,h]=onekpoiss(10);
scatter3(m,c,h)

for u=0.5:0.1:10
    [m,c,h]=onekpoiss(u,1);
    scatter3(m,c,h,5,'b')
end

function [m,c,h]=onekpoiss(u,n)
    if nargin<2
        n=1000;
    end
    c=zeros(n,1);
    m=zeros(n,1);
    h=zeros(n,1);
    for k=1:n
        a=poissrnd(u,[500,1]);
        m(k)=mean(a);
        c(k)=std(a)./mean(a);
        h(k)=sum(a==0);
    end
end

