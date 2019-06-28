function [We_cell] = simmaster_nodata(n, m, d, repprop, X_noise, S_noise, nsim)

%clc;clear all;

% Dataset Generation
%nsim = 100;

X_cell = cell(nsim,1);
W_cell = cell(nsim,1);

We_cell = cell(nsim,1);
MSE_L = zeros(nsim,1);
MSE_W = zeros(nsim,1);
MSE_W_nz = zeros(nsim,1);
MSE_W_nzzscore = zeros(nsim,1);
MSE_C = zeros(nsim,1);


We_cell_cor = cell(nsim,1);
MSE_L_cor = zeros(nsim,1);
MSE_W_cor = zeros(nsim,1);
MSE_W_nz_cor = zeros(nsim,1);
MSE_W_nzzscore_cor = zeros(nsim,1);
MSE_C_cor = zeros(nsim,1);

disp('Start running Algorithm');
for j = 1:nsim

      b = m/d; % Probes-Genes ratio
      X = zeros(n,m);
      W = zeros(m,d);
            for i = 1:d
                ind = (1+(i-1)*b):(i*b);
                if rand > repprop || i==1
                    sig = randn(n,1);
                end
                X(:,ind) = sig*ones(1,b) + X_noise*randn(n,b);
                W(ind,i) = rand(b,1);
            end
        X_cell{j} = X;
        W_cell{j} = W;
        L = X*W;
        S_true = corr(L);
        S_input = S_true + S_noise * randn(d,d); 
     
     
    % Run estW
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
        
        %run estW_cor
        We_cor = estW_cor(X,S_input,W);
        We_cell_cor{j} = We_cor;
        % Estimated L, and S
        Le_cor = X * We_cor;
        Se_cor = corr(Le_cor);

        % MSE: Average out all genes and samples

        MSE_L(j) = mean((L-Le_cor).^2,"All"); 
        MSE_W(j) = mean((W-We_cor).^2,"All"); 
        MSE_W_nz_cor(j) = mean((W(W~=0) - We_cor(W~=0)).^2);
        MSE_W_nzzscore_cor(j) = mean((zscore(W(W~=0)) - zscore(We_cor(W~=0))).^2);
        MSE_C(j) = sum((S_true-Se_cor).^2 .* ~eye(d),"All")/(d*(d-1)); 
        
    end
disp('Finish running algorithm');    
fname = sprintf('results/simres_all_n%d_m%d_d%d_repprop%.1f_Snoise%.1f_Xnoise%d_nsim%d.mat', n, m, d, repprop, S_noise, X_noise,nsim);
%'X_cell', 'W_cell', 'We_cell', 'MSE_L','MSE_C','MSE_W','MSE_W_nz', 'MSE_W_nzzscore','-v7.3'
save(fname);
disp('Success, Results saved!');
end