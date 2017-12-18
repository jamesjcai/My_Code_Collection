%--------------------------------------------------------------------------
% plotResults plots the results of the MIDER algorithm
% 
% Inputs: 
%   Output    = Structure that contains MIDER results
%   ntotal    = number of variables
%   variables = n-vector of strings with the names of the variables
%   options   = Structure that contains MIDER options
% 
% Outputs: one plot of the distances between variables and their links, and
% ntotal plots of the time-lagged mutual information between variables
%
%	Dependencies: plotResults.m calls projectVars.m, arrow.m, and requires
%	the Statistics toolbox (in Matlab) or the Statistics package (in Octave)
%   plotResults.m is called by runMIDER.m
%
% Written by Alejandro Fernández Villaverde (afvillaverde@iim.csic.es)
% Created: 22/03/2013, last modified: 29/01/2015
%--------------------------------------------------------------------------
% Copyright (C) 2014  Alejandro Fernández Villaverde
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------


function plotResults(Output,ntotal,variables,options)

% Pairwise distances:
tic
scaled_centered_dist = zscore(Output.dist);
dissimilarities      = pdist(scaled_centered_dist,'euclidean');
if(min(dissimilarities)<0.5)
    dissimilarities = dissimilarities + mean(dissimilarities);
    fprintf(1,' Some variables are very close. The distances between variables will be increased to improve visualization.');
end

% Project distances to lower-dimensional space using Multi-dimensional scaling:
Y = projectVars(dissimilarities,ntotal);
fprintf(1,' Multi-dimensional scaling performed in %d seconds.',toc); 

% Plot mutual information matrices:
if options.plotMI == 1
    if ntotal < 50
        [gridvars,gridtimes] = meshgrid(1:ntotal, -options.taumax:options.taumax );
        miarray = zeros(ntotal,options.taumax+1);
        miarrayneg = zeros(ntotal,options.taumax+1); 
        for i=1:ntotal
            figure(i)
            miarray(:,:)    = Output.MIl(i,:,:);
            miarrayneg(:,:) = Output.MIl(:,i,:);
            miarraynegorder = fliplr(miarrayneg);
            miarraycomplete = [miarraynegorder(:,1:end-1), miarray];
            surf(gridtimes, gridvars, miarraycomplete'); 
            titletext = sprintf('Mutual Information between %s and the other variables, at different lags', variables{i});
            title(titletext)   
            xlabel('Time lags')   
            ylabel('Variables')
            set(gca,'YTick',1:ntotal)
            set(gca,'YTickLabel',variables)
            zlabel('Normalized (0-1) Mutual Information')
        end
        figure(ntotal +1) 
    else
        fprintf(1,'\n MI plots were not displayed, to prevent memory errors due to the large number of variables (%d)',ntotal);
        fprintf(1,'\n If you still want to plot them, comment lines 54 & 74-77 in plotResults.m');
    end
end

% Plot map. First, draw arrows using function 'arrow.m':
axis(axis)
hold on
con_array = Output.con_array;
for i=1:ntotal
    for j=1:ntotal
        if con_array (i,j) > 0 
            if Output.T(i,j) > 0 %Output.T(j,i)
                arrow([Y(i,1),Y(i,2)], [Y(j,1),Y(j,2)],...
                    'Width', abs(50*(con_array(i,j))^1.5),'Length',90,'BaseAngle',40,'TipAngle',30)  
            else
                if Output.T(j,i) > 0
                    arrow([Y(j,1),Y(j,2)], [Y(i,1),Y(i,2)],...
                    'Width', abs(50*(con_array(i,j))^1.5),'Length',90,'BaseAngle',40,'TipAngle',30)
                else
                    arrow([Y(j,1),Y(j,2)], [Y(i,1),Y(i,2)],...
                    'Width', abs(50*(con_array(i,j))^1.5),'Length',0)
                end
            end
        end
    end
end


% Second, plot data points and names:
plot(Y(:,1),Y(:,2),'ok','Markersize',2,'MarkerFaceColor','k')
text(Y(:,1)+0.05,Y(:,2),variables,'FontSize',14)
minaxis = floor(min(min(Y(:,1:2))));
maxaxis = ceil(max(max(Y(:,1:2))));
axis([minaxis maxaxis minaxis maxaxis]); 
axis equal

end