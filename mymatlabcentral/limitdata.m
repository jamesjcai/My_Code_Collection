function [A2,B2]=limitdata(A,B,lima,limb)
%[A2,B2]=limitdata(A,B,lima,limb)

if nargin<3, error('3 para pls'); end

if(0)
B2=B(A<lima);
A2=A(A<lima);

if nargin>3
    A2=A2(B2<limb);
    B2=B2(B2<limb);
end
end


B2=B(A>lima);
A2=A(A>lima);

if nargin>3
    B2=B2(A2<limb);
    A2=A2(A2<limb);    
end



