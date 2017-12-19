function [E, p] = logistestt(x, y, w)
% labels = logistictest(x,y,w)
% computes the discriminative function based on log. reg. model
% with parameters w, returns the predicted labels


 p = g(x*w);
 
 E = y-p;
 
function p = g(z)
 p = 1./(1+exp(-z));
