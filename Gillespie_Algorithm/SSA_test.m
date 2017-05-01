%DEFINING REACTION RATES
par.kr = 1/10;     %transcription rate
par.kp = 10;   %translation rate
par.gr = 1/150; %mRNA degradation rate
par.gp = 1/30;  %protein degradation rate
 
%DEFINING PROPENSITY FUNCTIONS AS A FUNCTION HANDLE
prop = @(x,par)([par.kr,...      %transcription, one mRNA molecule created
                 par.kp*x(1),... %translation, one protein created
                 par.gr*x(1),... %mRNA degradation, one mRNA molecule removed
                 par.gp*x(2)]);  %protein degradation, one protein molecule removed
 
%DEFINING INITIAL CONDITION, order [mRNA, Protein]
init = [0;1];
 
%DEFINING STOICHIOMETRIC MATRIX
%column corresponds to the reaction, row corresponds to the molecule
%order as in prop and init variables
stoch = [1 0 -1 0;...
         0 1 0 -1];
 
%DEFINING TIME MESH FOR THE OUTPUT TRAJECTORY
tmesh = linspace(0,1000,100) ;


%simulating
traj = SSA(tmesh, par,prop,stoch, init );
 
%plotting
figure(1)
plot(tmesh,traj(1,:))
xlabel('Time, min')
ylabel('#mRNAs')
 
figure(2)
plot(tmesh,traj(2,:),'r')
xlabel('Time, min')
ylabel('#Proteins')
