addpath('fileexchange_14888-mutual-information-computation\');
a = [1 2 1 2 1]';
b = [2 1 2 1 1]';
c = [2 1 2 2 1]';
mutualinfo(a,b)
%rmpath('fileexchange_14888-mutual-information-computation\');

% WARNING!!! WARNING!!!! WARNING!!!
% This code assumes that you are working with discrete (integer) variables. If you are not, the code rounds your variables, turning them into integers without warning you. This will cause major errors if, for example, your variables are not integers and are less than 0.5. It will cause less extreme errors, but errors nonetheless, if your variables are larger, but not integers.
% For example the following code returns mutual information of zero! Which is obviously false.
x=[0 1 1 1 0 0 1 0 1 0 1 1 1 1 0 0 0 1 1]'*0.1; 
y=x; 
mutualinfo(x,y)

%%
addpath('MIToolbox-master\matlab\');
which mi
mi(a,b)
%rmpath('MIToolbox-master\matlab\');
mi(x,y)


%%
KSG_estimator_jc_parfor(x,y,5)
