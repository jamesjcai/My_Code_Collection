% ============================================================
% MDR DEMO: Detecting Epistasis via Multifactor Dimensionality Reduction
% Example: Pure XOR interaction (no marginal effects)
% ============================================================

clear; clc; rng(1)

%% ------------------------------------------------------------
% 1. Simulate epistatic genotype–phenotype data
% ------------------------------------------------------------
n = 400;

% Two SNPs coded as 0/1
SNP1 = randi([0 1], n, 1);
SNP2 = randi([0 1], n, 1);

% Pure epistasis: XOR interaction
Y = xor(SNP1, SNP2);   % 1 = case, 0 = control

fprintf('Marginal effects (should be ~0):\n')
fprintf('Corr(SNP1, Y) = %.3f\n', corr(SNP1, Y))
fprintf('Corr(SNP2, Y) = %.3f\n\n', corr(SNP2, Y))

%% ------------------------------------------------------------
% 2. Construct MDR contingency table
% ------------------------------------------------------------
case_count    = zeros(2,2);
control_count = zeros(2,2);

for i = 1:n
    if Y(i) == 1
        case_count(SNP1(i)+1, SNP2(i)+1) = ...
            case_count(SNP1(i)+1, SNP2(i)+1) + 1;
    else
        control_count(SNP1(i)+1, SNP2(i)+1) = ...
            control_count(SNP1(i)+1, SNP2(i)+1) + 1;
    end
end

disp('Case counts (rows=SNP1, cols=SNP2):')
disp(case_count)

disp('Control counts (rows=SNP1, cols=SNP2):')
disp(control_count)

%% ------------------------------------------------------------
% 3. MDR risk labeling (dimension reduction)
% ------------------------------------------------------------
threshold = 1;   % balanced case-control
risk_ratio = case_count ./ max(control_count, 1);
high_risk = risk_ratio > threshold;

disp('High-risk genotype combinations (1 = high risk):')
disp(high_risk)

%% ------------------------------------------------------------
% 4. Collapse multilocus genotypes → single MDR feature
% ------------------------------------------------------------
MDR_feature = zeros(n,1);
for i = 1:n
    MDR_feature(i) = high_risk(SNP1(i)+1, SNP2(i)+1);
end

%% ------------------------------------------------------------
% 5. Evaluate MDR prediction accuracy
% ------------------------------------------------------------
mdr_accuracy = mean(MDR_feature == Y);
fprintf('MDR accuracy = %.3f\n\n', mdr_accuracy)

%% ------------------------------------------------------------
% 6. Compare with logistic regression (no interaction term)
% ------------------------------------------------------------
X = [SNP1 SNP2];
b = glmfit(X, Y, 'binomial');
p = glmval(b, X, 'logit') > 0.5;

logit_accuracy = mean(p == Y);
fprintf('Logistic regression accuracy (no interaction) = %.3f\n\n', ...
        logit_accuracy)

%% ------------------------------------------------------------
% 7. Interpretation
% ------------------------------------------------------------
fprintf('Interpretation:\n')
fprintf(['- Individual SNPs show no marginal association.\n' ...
         '- MDR detects the joint XOR pattern by labeling\n' ...
         '  multilocus genotypes as high/low risk.\n' ...
         '- Logistic regression without interaction fails.\n'])

%% https://chatgpt.com/share/697ecab5-ada0-8005-bbb1-3d7d28d9db20