function i_scatter_withinfo(x,y,infotxt,X,Y)
if nargin<5, Y=[]; end
if nargin<4, X=[]; end
scatter(x,y,[],log(y)-log(x));
dt = datacursormode;
if isempty(X)
    dt.UpdateFcn = {@i_myupdatefcn1,infotxt};
else
    dt.UpdateFcn = {@i_myupdatefcn3,infotxt,X,Y};
end