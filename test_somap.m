
net = selforgmap([1 4]);
%view(net)
net = train(net,x);
%nntraintool
%plotsompos(net,x);
y = net(x);
cluster_indices = vec2ind(y);
%plotsomhits(net,x);
sc_scatter(X,g,s,cluster_indices)

% net = selforgmap([10 10]);
% [net,tr] = train(net,inputs); % Classify input data
% op_som = vec2ind(net(inputs))';


