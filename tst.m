clc; clear all;
%% Test Case: Data generation
n = 500;
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

L = X*W;
S = corr(L);
nfolds = 2;
alphas = [0:0.01:0.1,0.2,0.5,1];

mean_criterion_byalpha = zeros(length(alphas),1);
% Cross-validation Partition
CVO = cvpartition(n,'k',nfolds);
criterion_byfold = zeros(CVO.NumTestSets,1);

for k = 1:length(alphas)
    for i = 1:CVO.NumTestSets
        trIdx = CVO.training(i);
        teIdx = CVO.test(i);
        [Wl,criterion] = estW_smooth(X(trIdx,:),X(teIdx,:),S,W,alphas(k));
        criterion_byfold(i) = criterion;  
    end
    mean_criterion_byalpha(k) = mean(criterion_byfold);
    
end
[max_criterion, Idx_max_criteirion] = max(mean_criterion_byalpha);
best_alpha = alphas(Idx_max_criteirion);
figure;plot(alphas, mean_criterion_byalpha)

