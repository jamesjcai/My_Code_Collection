function [i]=denserank(a)



%a=[0 9 12 7; 3 0 8 2; 8 20 0 3; 25 2 18 0];
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/163003

needtransp=false;
if size(a,2)==1, a=a'; needtransp=true; end
[s i]=sort(a,2);
[J I]=ndgrid(1:size(i,1),1:size(i,2));
i(sub2ind(size(i),J,i))=I;
if needtransp, i=i'; end


%There exist many ranking functions, check here:
%http://en.wikipedia.org/wiki/Ranking.
%
%The most important are:
%Fractional Ranking (1 2.5 2.5 4)
%Dense Ranking (1 2 2 3)
%Standard Competition Ranking (1 2 2 4)
%Modified Competition Ranking (1 3 3 4)
%Ordinal Ranking (1 2 3 4) OR (1 3 2 4)
%****************************