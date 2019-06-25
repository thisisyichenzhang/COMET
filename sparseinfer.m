%% Test Case: Data generation
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
    X(:,ind) = sig*ones(1,b) + 1*randn(n,b);
    W(ind,i) = rand(b,1);
end

L = X*W;
S = corr(L);
[n,m] = size(X);
d = size(S,1);
% Initialization of W.
Wo = randn(m,d);
thr = 0.3;
Wl = estWnew(X,S,Wo,thr);






figure; 
subplot(1,3,1); imagesc(W)
subplot(1,3,2); imagesc(Wo)
subplot(1,3,3); imagesc(Wl)

L = X*W;
S = corr(L);
% initial
Lo = X * Wo;
So = corr(Lo);

% Estimated L, and S
Ll = X * Wl;
Sl = corr(Ll);

figure;
subplot(1,3,1);
imagesc(S.*~eye(d));
title("gene-wise corr (L)");
subplot(1,3,2);
imagesc(Sl.*~eye(d));
title("Gene wise corr (L_{l})");
subplot(1,3,3);
imagesc(So.*~eye(d));
title("Gene wise corr (L_{o})");


%ythard = wthresh(y,'h',thr);
%ytsoft = wthresh(y,'s',thr);
