function txt = i_myupdatefcn3(~,event_obj,g,X)
% Customizes text of data tips
% pos = event_obj.Position;
idx = event_obj.DataIndex;
% i_plotsiglegene(idx,g);
txt = {g(idx)};
figure;
plot(X(idx,:));
title(txt)
end