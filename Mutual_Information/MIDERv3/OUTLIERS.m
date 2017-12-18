%--------------------------------------------------------------------------
% OUTLIERS: outlier detection
%    
%   [ Xmd,OL ] = OUTLIERS( Xol ) takes as input a dataset Xol, which may 
%   contain outliers, and returns as output a modified dataset Xmd, where
%   the outliers have been removed (i.e. every outlier is replaced by a
%   NaN), and an array OL, which lists the locations of the outliers. 
%   Note that, after running OUTLIERS to detect the outlier values in a 
%   dataset, it is advisable to run the companion script TSR on the 
%   resulting dataset Xmd, which results in a new, corrected dataset with  
%   imputed values instead of the missing ones.
% 
%   Xol     = p*n array of input data (p data points; n variables)
%   
%   Xmd     = p*n array of output data (p data points; n variables)
%   OL      = positions (row,column) of the extreme outliers detected
%
%	Dependencies: When used as part of the MIDER toolbox, OUTLIERS.m is 
%   called by runMIDER.m
% 
%   Written by Abel Folch Fortuny (abfolfor@upvnet.upv.es)
%   Created: 01/12/2014, last modified: 27/01/2015
%--------------------------------------------------------------------------
% Copyright (C) 2015  Abel Folch Fortuny
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

function [ Xmd,OL ] = OUTLIERS( Xol )

[m,n]=size(Xol);

x_as=(Xol-repmat(mean(Xol),m,1))./repmat(std(Xol),m,1);

[a,b,c]=svd(cov(x_as));
diagb=diag(b);
ncomp=sum(diagb>=1);

[uu,ss,vv]=svd(x_as);
xrec=uu(:,1:ncomp)*ss(1:ncomp,1:ncomp)*(vv(:,1:ncomp))';
diag(ss);

E=x_as-xrec;
EE=E*E';
EE=diag(EE);

perclim=0.95;

numb=max(1,floor(perclim*m));
rest=m-numb;

fil=0;
for round=1:1000
    
    ord=randperm(m); 
    subset=EE(ord(1:numb));
    fil=fil+1;
    [d I]=sort(subset,'descend');
    LIM(fil)=subset(I(rest+1))+0.001;

end   
LIM=median(LIM);

% OUTLIERS
OL=[];
[d I]=sort(EE,'descend');
for i=1:length(EE)
    if d(i)>LIM;
        OL=[OL; I(i)];
    end
end
L=length(OL);

% 2 times above LIMIT
OL_2times=[];
i=1;
stop=0;
while stop<1 && i<=length(OL)
   
    if EE(OL(i))>2*LIM
       OL_2times=OL(1:i); 
       i=i+1;
    else
       stop=1;
    end
    
end

% 10 times above LAST
OLnew=[];
EE_new=EE(OL(1:L))-ones(L,1)*LIM;
LAST=EE_new(end);

i=L;
stop=0;
while stop<1 && i>=2
   
    i=i-1;
    if EE_new(i)<10*LAST
       OLnew=OL(1:(i-1));  
    else
        stop=1;
    end
    
end

OL=unique([OLnew;OL_2times]);
L=length(OL);

if L>0
    for l=1:L

        E2_OL=E(OL(l),:).^2;
        [d I]=sort(E2_OL,'descend');
        OL(l,2)=I(1);

    end
end

Xmd=Xol;

if size(OL,1)>0
    Xmd(OL(:,1),OL(:,2))=nan;
end
