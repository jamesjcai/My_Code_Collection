function i_scatter_withinfo(x,y,infotxt,X)
if nargin<4, X=[]; end
scatter(x,y);
dt = datacursormode;
if isempty(X)
    dt.UpdateFcn = {@i_myupdatefcn1,infotxt};
else
    dt.UpdateFcn = {@i_myupdatefcn3,infotxt,X};
end