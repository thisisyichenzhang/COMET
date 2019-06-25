
function [W_cell, We_cell,X_cell, MSE_L, MSE_W, MSE_W_nz, MSE_W_nzzscore, MSE_C] = simmaster(n, m, d, X_noise, S_noise, nsim)

%clc;clear all;

% Dataset Generation
%nsim = 100;
W_cell = cell(nsim,1);
We_cell = cell(nsim,1);
X_cell = cell(nsim, 1);
S_cell = cell(nsim, 1);
MSE_L = zeros(nsim,1);
MSE_W = zeros(nsim,1);
MSE_W_nz = zeros(nsim,1);
MSE_W_nzzscore = zeros(nsim,1);
MSE_C = zeros(nsim,1);

for j = 1:nsim
% d/n  = 100/100; 500/100; 500/1000
%m = 3000 and 15000 for d = 100 and 500 so that m/d ~= 30.
%n = 1000; % Number of subjects
%m = 100; % Number of Probes
%d = 20;  % Number of genes
%X_noise = 10 
b = m/d; % Probes-Genes ratio
X = zeros(n,m);

W = zeros(m,d);
    for i = 1:d
        ind = (1+(i-1)*b):(i*b);
        if rand > 0.25 || i==1
            sig = randn(n,1);
        end
        X(:,ind) = sig*ones(1,b) + X_noise*randn(n,b);
        W(ind,i) = rand(b,1);
    end
X_cell{j} = X;
W_cell{j} = W;

% Run
    L = X*W;
    S_true = corr(L);
    S_input = S_true + S_noise * randn(d,d); 
    We = estW(X,S_input,W);
    We_cell{j} = We;
    % Estimated L, and S
    Le = X * We;
    Se = corr(Le);

    % MSE: Average out all genes and samples

    MSE_L(j) = mean((L-Le).^2,"All"); 
    MSE_W(j) = mean((W-We).^2,"All"); 
    MSE_W_nz(j) = mean((W(W~=0) - We(W~=0)).^2);
    MSE_W_nzzscore(j) = mean((zscore(W(W~=0)) - zscore(We(W~=0))).^2);
    MSE_C(j) = sum((S_true-Se).^2 .* ~eye(d),"All")/(d*(d-1)); 
end
fname = sprintf('myfile_n%d_m%d_d%d_Snoise%.1f_Xnoise%d.mat', n, m, d, S_noise, X_noise);
save(fname);
end





