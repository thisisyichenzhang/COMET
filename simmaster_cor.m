function [We_cell] = simmaster_cor(n, m, d, repprop, X_noise, S_noise, nsim)

%clc;clear all;
dataname = sprintf('data/data_n%d_m%d_d%d_repprop%.1f_Snoise%.1f_Xnoise%d_nsim100.mat', n, m, d, repprop, S_noise, X_noise);
load(dataname);
disp('Load the data successfully');
% Dataset Generation
%nsim = 100;

We_cell = cell(nsim,1);
MSE_L = zeros(nsim,1);
MSE_W = zeros(nsim,1);
MSE_W_nz = zeros(nsim,1);
MSE_W_nzzscore = zeros(nsim,1);
MSE_C = zeros(nsim,1);

disp('Start running estW_cor');
for j = 1:nsim

     X = X_cell{j};
     W = W_cell{j};
     S_input = S_input_cell{j};
     S_true = S_true_cell{j};
    % Run estW_cor
        L = X*W;

        We = estW_cor(X,S_input,W);
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
disp('Finish running algorithm');    
fname = sprintf('results/simres_cor_n%d_m%d_d%d_repprop%.1f_Snoise%.1f_Xnoise%d_nsim%d.mat', n, m, d, repprop, S_noise, X_noise,nsim);
save(fname,'We_cell', 'MSE_L','MSE_C','MSE_W','MSE_W_nz', 'MSE_W_nzzscore','-v7.3');
disp('Success, Results saved!');
end