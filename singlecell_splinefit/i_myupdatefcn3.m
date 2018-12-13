function txt = i_myupdatefcn3(~,event_obj,g,X,Y)
if nargin<5
    Y=[];
end
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
figure;
stem(1:length(X(idx,:)), X(idx,:),'marker','none');
if ~isempty(Y)
     hold on
     stem(1+length(X(idx,:)):length(Y(idx,:))+length(X(idx,:)),...
          Y(idx,:),'marker','none');
end
title(txt)
end