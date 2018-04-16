function [X]=s_square_diff(y)
% n=length(y);
% X=zeros(n);
% 
% for i=1:n-1
%     for j=i+1:n
%         %X(i,j)=(y(i)-y(j)).^2;
%         X(i,j)=(y(i)-y(j))*(y(i)-y(j));
%     end
% end
X=(y-y').^2;