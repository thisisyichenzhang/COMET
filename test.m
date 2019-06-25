clear all; clc;
% Test Case: 

% Data generation
n = 1000;
m = 100;
d = 20;
b = m/d;
X = zeros(n,m);
W = zeros(m,d);
for i = 1:d
    ind = (1+(i-1)*b):(i*b);
    if rand > 0.25 || i==1
        sig = randn(n,1);
    end
    X(:,ind) = sig*ones(1,b) + 0.1*randn(n,b);
    W(ind,i) = rand(b,1);
end

% How does true W looks like?
figure; imagesc(W)
%figure; imagesc(corr(X).*~eye(size(X,2)))
figure; imagesc(corr(X))
K = X * W;
%figure; imagesc(corr(K).*~eye(size(K,2)))
figure; imagesc(corr(K))


% Run the algorihtm 
L = X*W;
S = corr(L);
Wl = estW_cor(X,S,W);

% How does estiameted W look like?
figure; 
subplot(1,2,1); imagesc(zscore(Wl,[],1)); title("Estimated  (Z-scored)");
subplot(1,2,2); imagesc(zscore(W,[],1)); title("True W (Z-scored)");

figure; plot(zscore(W(W~=0)));
hold on; plot(zscore(Wl(W~=0)));

Khat = X * Wl;
figure; 
subplot(1,2,1); imagesc(corr(Khat)); title("Estimated correlations for gene expression");
subplot(1,2,2); imagesc(corr(K)); title("True correlations for gene expression");
% How does the estimated correlations across genes compared to the given 
