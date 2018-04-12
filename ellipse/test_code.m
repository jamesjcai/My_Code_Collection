load testdata datax xy

fig=figure;
plot(datax(:,1),datax(:,2),'o')
hold on
error_ellipse(cov(datax),mean(datax),'conf',0.95);
plot(xy(1),xy(2),'*r');

plot(0.8182,0.066,'*');

dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn',{@myupdatefcn,datax})

%%

% this should return FALSE
i_95pct_error_ellipse(datax,xy,0.95)

% this should return TRUE
i_95pct_error_ellipse(datax,mean(datax),0.95)

% this should return TRUE
i_95pct_error_ellipse(datax,[0.8182 0.066],0.95)


function txt = myupdatefcn(~,event_obj,datax)
% Customizes text of data tips
pos = get(event_obj,'Position');
%I = get(event_obj, 'DataIndex');
yes=i_95pct_error_ellipse(datax,pos,0.95);
if yes
    t='inside';
else
    t='outside';
end
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       t};
end
