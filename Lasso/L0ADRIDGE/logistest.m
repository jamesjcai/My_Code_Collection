function [labels, p] = logistest(x, w)
% labels = logistictest(x,y,w)
% computes the discriminative function based on log. reg. model
% with parameters w, returns the predicted labels


 p = g(x*w);
 
 labels = (p>0.5);
 
function p = g(z)
 p = 1./(1+exp(-z));
