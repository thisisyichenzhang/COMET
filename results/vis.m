clc; clear all;
load('myfile_n100_m3000_d100_Snoise0.0_Xnoise10.mat')
x1 = MSE_C; x1_MSE_L = MSE_L; x1_MSE_W = MSE_W; x1_MSE_W_nz = MSE_W_nz; x1_MSE_W_nzscore = MSE_W_nzzscore;
load('myfile_n100_m3000_d100_Snoise0.5_Xnoise10.mat')
x2 = MSE_C; x2_MSE_L = MSE_L; x2_MSE_W = MSE_W; x2_MSE_W_nz = MSE_W_nz; x2_MSE_W_nzscore = MSE_W_nzzscore;
load('myfile_n100_m3000_d100_Snoise1.0_Xnoise10.mat')
x3 = MSE_C; x3_MSE_L = MSE_L; x3_MSE_W = MSE_W; x3_MSE_W_nz = MSE_W_nz; x3_MSE_W_nzscore = MSE_W_nzzscore;
load('myfile_n100_m3000_d100_Snoise10.0_Xnoise10.mat')
x4 = MSE_C; x4_MSE_L = MSE_L; x4_MSE_W = MSE_W; x4_MSE_W_nz = MSE_W_nz; x4_MSE_W_nzscore = MSE_W_nzzscore;

% MSE_C
f1 = [x1;x2;x3;x4];
g = [0 .* ones(size(x1)); 0.5.*ones(size(x2)); ones(size(x3)); 10.*ones(size(x4))];
figure
boxplot(f1,g)
xlabel('Noise level (SD) of input matrix K ')
ylabel('MSE')
title('MSE for correlations estiamtes at different noise levels of K')
savefig("MSE_C");

% MSE_W
f2 = [x1_MSE_W;x2_MSE_W;x3_MSE_W;x4_MSE_W];
figure
boxplot(f2,g)
xlabel('Noise level (SD) of input matrix K ')
ylabel('MSE of W')
title('MSE for W estiamtes at different noise levels of K')
savefig("MSE_W");

% MSE_W_nz
f3 = [x1_MSE_W_nz;x2_MSE_W_nz;x3_MSE_W_nz;x4_MSE_W_nz];
figure
boxplot(f3,g)
xlabel('Noise level (SD) of input matrix K ')
ylabel('MSE of W')
title('MSE for W (nonzero) estiamtes at different noise levels of K')
savefig("MSE_W_nz");

% MSE_W_nz
f4 = [x1_MSE_W_nzscore;x2_MSE_W_nzscore;x3_MSE_W_nzscore;x4_MSE_W_nzscore];
figure
boxplot(f4,g)
xlabel('Noise level (SD) of input matrix K ')
ylabel('MSE of W')
title('MSE for non-zero W (zscore) at different noise levels of K')
savefig("MSE_W_nzzscore");

% MSE_L
f5 = [x1_MSE_L;x2_MSE_L;x3_MSE_L;x4_MSE_L];
figure
boxplot(f5,g)
xlabel('Noise level (SD) of input matrix K ')
ylabel('MSE of L')
title('MSE for L estimate (gene expression matrix) at different noise levels of K')
savefig("MSE_L");
