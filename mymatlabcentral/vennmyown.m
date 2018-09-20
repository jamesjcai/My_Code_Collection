function vennmyown(k1,k2,k3)

if islogical(k1)
    a7=sum(k1&k2&k3);
    a2=sum(and(k1,k2))-a7;
    a6=sum(and(k1,k3))-a7;
    a1=sum(k1)-a2-a6-a7;
    a4=sum(and(k2,k3))-a7;
    a3=sum(k2)-a2-a4-a7;
    a5=sum(k3)-a4-a6-a7;
else
    a7=length(intersect(k1,intersect(k2,k3)));
    a2=length(intersect(k1,k2))-a7;
    a6=length(intersect(k1,k3))-a7;
    a1=length(k1)-a2-a6-a7;
    a4=length(intersect(k2,k3))-a7;
    a3=length(k2)-a2-a4-a7;
    a5=length(k3)-a4-a6-a7;
end
%figure;
[ a1 a2 a3 a4 a5 a6 a7 ]

venn([ a1 a2 a3 a4 a5 a6 a7 ])
%vennX( [ a1 a2 a3 a4 a5 a6 a7 ], .1);
axis equal
