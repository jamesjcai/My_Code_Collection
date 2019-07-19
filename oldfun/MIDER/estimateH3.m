%--------------------------------------------------------------------------
% estimateH3: obtain joint entropy of three variables
% 
%   [H3,fracn] = estimateH3(x,y,z,pb,q) calculates the joint entropy
%   'H3' of three variables 'x', 'y', and 'z', using the type of entropy 
%   specified by parameter 'q', and estimating the entropy with 'pb' data 
%   points per bin. 
% 
%   x       = data vector of the 1st variable
%   y       = data vector of the 2nd variable 
%   z       = data vector of the 3rd variable 
%   pb      = estimated number of points per bin
%   q       = entropic parameter; if q>1 it is Tsallis entropy
%   fracn   = fraction of joint prob cells that contain >= 5 points.
%   H3      = joint entropy of (x,y,z)
% 
%   Dependencies: estimateH3.m is called by mider.m
%
%	Written by Alejandro Fernández Villaverde (afvillaverde@iim.csic.es)
%	Created: 07/10/2011, last modified: 24/11/2014
%--------------------------------------------------------------------------
% Copyright (C) 2014  Alejandro Fernández Villaverde
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

function [H3,fracn] = estimateH3(x,y,z,pb,q)

szx     = size(x); if szx(1)>szx(2), x=x'; end  
szy     = size(y); if szy(1)>szy(2), y=y'; end
szz     = size(z); if szz(1)>szz(2), z=z'; end
leng    = length(x);
binsize = ceil(sqrt(pb*leng));   % Estimated bin size for ~pb points/bin
b       = floor(leng/binsize);   % Number of bins
binsize = floor(leng/b);         % Integer number of elements per bin

% Initialize matrices to zero.
his  = zeros(b,b);
his3 = zeros(b,b,b);
x1   = zeros(b,1);
x2   = zeros(b,1);
x3   = zeros(b,1);

% Sort values of X, Y, and Z (descending order)
[sx,ix] = sort(x,2,'descend');    
[sy,iy] = sort(y,2,'descend');   
[sz,iz] = sort(z,2,'descend');   

% Replace elements of X, Y, and Z by their rank order
linsp  = linspace(1,leng,leng);
sx(ix) = linsp;
sy(iy) = linsp;
sz(iz) = linsp;

% Put sorted elements into bins
sx = ceil(sx/binsize);
sy = ceil(sy/binsize);
sz = ceil(sz/binsize);

% Put largest values in highest bin
for p=1:leng
   if sx(p) > b
        sx(p)=b;
   end
   if sy(p) > b
        sy(p)= b;
   end
   if sz(p) > b
        sz(p)= b;
   end
end

% Populate histograms
for i=1:leng,
    x = sx(i);
    y = sy(i);
    z = sz(i);
    his(y,x)    = (his(y,x)+1);
    his3(z,y,x) = (his3(z,y,x)+1);
    x1(x) = (x1(x)+1);
    x2(y) = (x2(y)+1);
    x3(z) = (x3(z)+1);
end
                    
% Continue with calculation of H3
his  = his/leng;
his3 = his3/leng;
H3   = 0;
if q>1
    for j = 1:b,
        for k = 1:b,
            for l = 1:b,
                if(his3(j,k,l)>0)
                    an = - his3(j,k,l)*(his3(j,k,l)^(1-q)-1)/(1-q);
                    H3 = H3 + an;
                end
            end
        end
    end
else
    for j = 1:b,
        for k = 1:b,
            for l = 1:b,
                if(his3(j,k,l)>0)
                    an = - his3(j,k,l)*log2(his3(j,k,l));   
                    H3 = H3 + an;
                end
            end
        end
    end
end

% Calculate fracn
[a,b]   = size(his);
his     = his*leng;
m =1 ; 
for k=1:a;
    for j=1:b;
        if his(k,j)>0,
            hnz(m)=his(k,j);
            m=m+1;
        end;
    end;
end
lenz = length(hnz); 
num5 = 0; for k=1:lenz;if hnz(k)>=5;num5=num5+1;end;end
fracn = num5/lenz;
