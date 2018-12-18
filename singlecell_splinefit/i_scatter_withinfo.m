function i_scatter_withinfo(d0,d1,infotxt,X0,X1)
if nargin<5, X1=[]; end
if nargin<4, X0=[]; end
scatter(d0,d1,[],log(d1)-log(d0));
dt = datacursormode;
if isempty(X0)
    dt.UpdateFcn = {@i_myupdatefcn1,infotxt};
else
    dt.UpdateFcn = {@i_myupdatefcn3,infotxt,X0,X1};
end