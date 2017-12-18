%--------------------------------------------------------------------------
% TSR: Trimmed Scores Regression
%    
%   [X]     = TSR(Xmd) takes as input a dataset Xmd, which may be 
%   incomplete (i.e. it may have missing values), and returns as output a 
%   complete dataset X, where the missing entries have been replaced with
%   imputed values computed with the Trimmed Scores Regression method.
% 
%   Xmd     = p*n array of input data (p data points; n variables)
%   
%   X       = p*n array of output data (p data points; n variables)
%
%	Dependencies: When used as part of the MIDER toolbox, TSR.m is called 
%   by runMIDER.m
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

function [X]=TSR(Xmd)

[n,p]=size(Xmd);

ve=0;
ncomp=0;
if exist('octave_config_info','builtin') == 0
    B=corr(Xmd,'rows','pairwise');  % Matlab command
else
    B=corrcoef(Xmd);                % Octave command
end
B(isnan(B)==1)=0;
[u,s,v]=svd(B);
diags=diag(s);
while ve<0.9
    ncomp=ncomp+1;
    ve=ve+diags(ncomp)/sum(diags);
end

for i=n:-1:1,
  r=~isnan(Xmd(i,:));
  pat(i).O=find(r==1); 
  pat(i).M=find(r==0); 
  pat(i).nO=size(pat(i).O,2); 
  pat(i).nM=size(pat(i).M,2); 
end

X=Xmd;
mis=isnan(Xmd);
[r c]=find(isnan(Xmd));
X(mis)=0;
meanc=sum(X)./(n-sum(mis));
for k=1:length(r),
  X(r(k),c(k))=meanc(c(k));
end

maxiter=5000;
conv=1.0e-10;
diff=100;

It=0;
while It<maxiter && diff>conv,
  It=It+1;
  Xmis=X(mis);
  mX=mean(X);
  S=cov(X);
  Xc=X-ones(n,1)*mX;
  if n>p, [U D V]=svd(Xc,0); else [V D U]=svd(Xc',0); end

  V=V(:,1:ncomp);
  for i=1:n,              % for each row
    if pat(i).nM>0,       % if there are missing values
      L=V(pat(i).O,1:min(ncomp,pat(i).nO));
      S11=S(pat(i).O,pat(i).O);
      S21=S(pat(i).M,pat(i).O);
      z1=Xc(i,pat(i).O)';
      z2=S21*L*pinv(L'*S11*L)*L'*z1;  % use pseudoinverse
      Xc(i,pat(i).M)=z2';
    end
  end
  X=Xc+ones(n,1)*mX;
  d=(X(mis)-Xmis).^2;
  diff=mean(d);
end

