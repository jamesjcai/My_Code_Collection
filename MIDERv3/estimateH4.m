%--------------------------------------------------------------------------
% estimateH4: obtain joint entropy of four variables
% 
%   [H4,fracn] = estimateH4(w,x,y,z,pb,q) calculates the joint entropy
%   'H4' of four variables 'w', 'x', 'y', and 'z', using the type of entropy 
%   specified by parameter 'q', and estimating the entropy with 'pb' data 
%   points per bin. 
% 
%   w       = data vector of the 1st variable
%   x       = data vector of the 2nd variable
%   y       = data vector of the 3rd variable 
%   z       = data vector of the 4th variable 
%   pb      = estimated number of points per bin
%   q       = entropic parameter; if q>1 it is Tsallis entropy
%   fracn   = fraction of joint prob cells that contain >= 5 points.
%   H4      = joint entropy of (w,x,y,z)
% 
%   Dependencies: estimateH4.m is called by mider.m
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

function [H4,fracn] = estimateH4(w,x,y,z,pb,q)

szw     = size(w); if szw(1)>szw(2), w=w'; end
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
his4 = zeros(b,b,b,b);
x0   = zeros(b,1);
x1   = zeros(b,1);
x2   = zeros(b,1);
x3   = zeros(b,1);

% Sort values of W, X, Y, and Z (descending order)
[sw,iw] = sort(w,2,'descend');   
[sx,ix] = sort(x,2,'descend');     
[sy,iy] = sort(y,2,'descend');   
[sz,iz] = sort(z,2,'descend');   

% Replace elements of W, X, Y and Z by their rank order
linsp  = linspace(1,leng,leng);
sw(iw) = linsp;
sx(ix) = linsp;
sy(iy) = linsp;
sz(iz) = linsp;

% Put sorted elements into bins
sw = ceil(sw/binsize);
sx = ceil(sx/binsize);
sy = ceil(sy/binsize);
sz = ceil(sz/binsize);

%Put largest values in highest bin
for p=1:leng
   if sw(p) > b
       sw(p)=b;
   end
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

%Populate histograms
for i=1:leng,
    w = sw(i);
    x = sx(i);
    y = sy(i);
    z = sz(i);
    his(y,x)    = his(y,x)+1;
    his3(z,y,x) = his3(z,y,x)+1;
    his4(z,y,x,w) = his4(z,y,x,w)+1;
    x0(w) = x0(w)+1;
    x1(x) = x1(x)+1;
    x2(y) = x2(y)+1;
    x3(z) = x3(z)+1;
end

% Continue with calculation of H4
his  = his/leng;
his4 = his4/leng;
H4   = 0;
if q>1
    for j = 1:b,
        for k = 1:b,
            for l = 1:b,
                for i = 1:b,
                    if(his4(j,k,l,i)>0)
                        an = - his4(j,k,l,i)*(his4(j,k,l,i)^(1-q)-1)/(1-q);
                        H4 = H4 + an;
                    end
                end
            end
        end
    end
else
    for j = 1:b,
        for k = 1:b,
            for l = 1:b,
                for i = 1:b,
                    if(his4(j,k,l,i)>0)
                        an = - his4(j,k,l,i)*log2(his4(j,k,l,i));   % /(x1(k).*x2(j).*x3(l)));
                        H4 = H4 + an; 
                    end
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
