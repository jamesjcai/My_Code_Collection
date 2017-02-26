function [yo,no,xo] = datapooling(y,x)
%
%   N = HIST(Y,M), where M is a scalar, uses M bins.
%
%   N = HIST(Y,X), where X is a vector, returns the distribution of Y
%   among bins with centers specified by X. The first bin includes
%   data between -inf and the first center and the last bin
%   includes data between the last bin and inf. Note: Use HISTC if
%   it is more natural to specify bin edges instead. 
%
%   [N,X] = HIST(...) also returns the position of the bin centers in X.
%
%   HIST(...) without output arguments produces a histogram bar plot of
%   the results. The bar edges on the first and last bins may extend to
%   cover the min and max of the data unless a matrix of data is supplied.
%
%   HIST(AX,...) plots into AX instead of GCA.
%
%   Class support for inputs Y, X: 
%      float: double, single
%
%   See also HISTC, MODE.

%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.20.4.9 $  $Date: 2005/12/12 23:24:07 $

% Parse possible Axes input


if nargin == 1, x = 10; end
if min(size(y))==1, y = y(:); end
if ~ishistnumeric(x) || ~ishistnumeric(y)
    error('MATLAB:hist:InvalidInput', 'Input arguments must be numeric.')
end

if isempty(y),
    if length(x) == 1,
       x = 1:double(x);
    end
    nn = zeros(size(x)); % No elements to count
    %  Set miny, maxy for call to bar below.
    miny = [];
    maxy = [];
else
    %  Ignore NaN when computing miny and maxy.
    ind = ~isnan(y);
    miny = min(y(ind));
    maxy = max(y(ind));
    %  miny, maxy are empty only if all entries in y are NaNs.  In this case,
    %  max and min would return NaN, thus we set miny and maxy accordingly.
    if (isempty(miny))
      miny = NaN;
      maxy = NaN;
    end
    if length(x) == 1
    	  if miny == maxy,
    		  miny = miny - floor(x/2) - 0.5; 
    		  maxy = maxy + ceil(x/2) - 0.5;
     	  end
        binwidth = (maxy - miny) ./ x;
        xx = miny + binwidth*(0:x);
        xx(length(xx)) = maxy;
        x = xx(1:length(xx)-1) + binwidth/2;
    end
    % Shift bins so the interval is ( ] instead of [ ).
    xx = full(real(xx)); y = full(real(y)); % For compatibility
    bins = xx + eps(xx);    
    [nn,yo] = histc(y,[-inf bins],1);

    yo(find(yo==1))=2;
    yo(find(yo==length(nn)))=length(nn)-1;
    yo=yo-1;
    
    % Combine first bin with 2nd bin and last bin with next to last bin
    nn(2,:) = nn(2,:)+nn(1,:);
    nn(end-1,:) = nn(end-1,:)+nn(end,:);
    nn = nn(2:end-1,:);
    
end

if nargout == 0
    bar(x,nn,[miny maxy],'hist');
else
  if min(size(y))==1, % Return row vectors if possible.
    no = nn';
    xo = x;
  else
    no = nn;
    xo = x';
  end
end

function a = ishistnumeric(b)
% for backward compatibility, logical is allowed in hist.m
a = isnumeric(b) || islogical(b);
