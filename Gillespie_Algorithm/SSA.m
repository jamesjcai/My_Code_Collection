function trajectory = SSA(tmesh, par,prop,stoch, init )
   %tmesh - time mesh on which solution should be returned
   %per - parameters of the pathway
   %prop - definition of propensity functions
   %stoch - stochiometric matrix
   %init - initial condition for the pathway
 
   t = 0;                                          %current time
   state = init(:);                                %variable with current system state
   trajectory = zeros(length(init),length(tmesh)); %preparing output trajectory
   trajectory(:,1) = init(:);                      %setting initial value as the first element in trajectory
   cindx = 2;                                      %current trajectory index
   N = length(tmesh);                              %number of time points
 
   while t<tmesh(end)
      Q = feval(prop,state,par);          %calculating propensities of the reactions
      Qs = sum(Q);                        %total propensity
      dt = -log(rand())/Qs;               %generating time to the next reaction
      R = sum(rand >= cumsum([0 Q])/Qs);  %selecting reaction
      state = state + stoch(:,R);         %updating state
      t = t + dt;                         %updating time
 
      %writing the output
      while cindx<=N && t>tmesh(cindx)
         trajectory(:,cindx) = state;
         cindx = cindx+1;
     end
   end
end
