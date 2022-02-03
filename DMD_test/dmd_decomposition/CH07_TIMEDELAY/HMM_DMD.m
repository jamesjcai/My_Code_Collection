clear all, close all, clc
% Weather state [S R C]: sun (S), rain (R), clouds (C)
% if sunny today: 0% sun tomorrow, 40% rain, 60% cloud
% rainy today: 25% sun tomorrow, 50% rain, 25% cloud
% cloudy today: 25% sun tomorrow, 25% rain, 50% cloud
% Transition Matrix T(i,j) = Prob(i|j), 
T = [0.00 0.40 0.60;
    0.25 0.50 0.25;
    0.25 0.25 0.50];
% x_{k+1} = x_k * T,     sum(T') = [1 1 1]

% Emissions E
% my dog is either happy (H) or grumpy (G)
% if it is sunny, dog is happy
% if it is rainy, dog is grumpy
% if it is cloudy, dog is %20 happy, dog is grumpy 80%
E = [1.0 0.0;
    0.0 1.0;
    0.2 0.8];

% Look at eigenvalues
[V,D] = eigs(T');  % transpose for left eigenvectors
V = V/sum(V(:,1)); % normalize probabilities = 1

%%  Compute a ten day forcast 
xToday = [0 1 0]; % day 0
xHist = xToday;   % keep track of all days
for i=1:10    % ten day forcast
    xTomorrow = xToday*T;
    xHist = [xHist; xTomorrow];
    xToday = xTomorrow;
end

figure
subplot(2,1,1)
plot(0:1:length(xHist)-1,xHist)
legend('Sunny','Rainy','Cloudy')
xlabel('# of Days from Today'), ylabel('Probability')
subplot(2,1,2)
plot(0:1:length(xHist)-1,xHist*E)
legend('Happy Dog','Grumpy Dog')
xlabel('# of Days from Today'), ylabel('Probability')

%% HMMGENERATE:  Generate states from T, E
N = 1000;
[seq,states] = hmmgenerate(N,T,E);
%     'Symbols',{'Happy','Grumpy'},...
%     'Statenames',{'Sunny';'Rainy';'Cloudy'})
indS = (states==1);
indR = (states==2);
indC = (states==3);
indH = (seq==1);
indG = (seq==2);

% HMMVITERBI: Probable states from Emissions, T, E
likelystates = hmmviterbi(seq, T, E);
percentCorrect = sum(likelystates==states)/N;

% HMMESTIMATE: Estimate T, E from Emissions, States
[T_est, E_est] = hmmestimate(seq, states);

% HMMDECODE: Posterior state probability from Emissions, T, E
pstates = hmmdecode(seq, T, E);

% HMMTRAIN: Estimate T, E from Emissions, T_guess, E_guess... algorithm isn't great for our data
T_guess = T + .001*randn(size(T));
E_guess = E + .001*randn(size(E));
[T_est2, E_est2] = hmmtrain(seq,T_guess,E_guess,'algorithm','viterbi');


%% FIGURES
figure
plot(states,'k-'), hold on
plot(likelystates,'rx')

 
 



%% Observability...
rank(obsv(T',E'))
E_bad = [.75 .25;
    .5 .5; 
    .5 .5];
rank(obsv(T',E_bad'))

[seq,states] = hmmgenerate(N,T,E_bad);%,...
likelystates = hmmviterbi(seq, T, E_bad);
[T_est, E_est] = hmmestimate(seq, states);

%% Connection with DMD... 
%    ... use ERA for original emmissions matrix E
E = eye(3);
N = 10000;
[seq,states] = hmmgenerate(N,T,E); 

indS = (states==1);
indR = (states==2);
indC = (states==3);

[T_est, E_est] = hmmestimate(seq, states);

X = zeros(3,N);
for i=1:N
    X(seq(i),i) = 1;
end

Y = X(:,2:end);
X = X(:,1:end-1);
[U,S,V] = svd(X,'econ');

Sinv = S^(-1);
Atilde = U(:,1:3)'*Y*V(:,1:3)*Sinv;
% step 3
[W,D] = eig(Atilde);
% step 4
Phi = Y*V(:,1:3)*Sinv*W;