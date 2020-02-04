load data/ROSMAP/exprByChr/expressionNormalizedChr18.mat;
load data/ROSMAP/methyByChr/methylationSNMnormChr18.mat;
load ('data/ROSMAP/mapping/methyToGeneByChr/methyToGeneChr18_1MB.mat');
n = 702; m = 5404; d = 250;

XtoY;

upperlimit = 5000;
lowerlimit = -5000;

W_mask = double( XtoY~=0 & XtoY<upperlimit & XtoY>lowerlimit);
%figure;
%subplot(121);stem(sum(W_mask)); title('Number of neighboring genes for each probe')
%subplot(122);stem(sum(W_mask,2)); title('Number of probes for each gene')

noise = 0;
W_new = rand(m,d) .* W_mask + noise * rand(m,d);
figure;
XtoY_aftermask = XtoY.* W_mask;
subplot(121);imagesc(XtoY_aftermask);
subplot(122);imagesc(W_new);

%%%%%%%%
n= 702;
m = 1832;
d = 198;
b= m/d;
repprop = 0.5;
nsim = 10;
seed = 123;
W_maskmat = load('W_mask.mat');
W_mask = W_maskmat.W_mask;
W_mask_active = W_mask(sum(W_mask~=0,2)~=0,sum(W_mask~=0,1)~=0);
rng(123);[X,W,L,W_init,L_input,X_input,S_input,paras]=simdata(n,m,d,0.6,0,0,0, W_mask_active);

figure;
subplot(221);imagesc(corr(X_input));title('Simulated co-methylation Network')
subplot(222);imagesc(corr(methy(:,sum(W_mask~=0,2)~=0)));title('Co-methylation Network - ROSMAP Chr18')
subplot(223);imagesc(corr(L_input));title('Simulated co-expression Network')
subplot(224);imagesc(corr(expr(:,sum(W_mask~=0,1)~=0)));title('Co-expression Network - ROSMAP Chr18')
save(gcf,"sim_ROSMAP_Chr18.jpg");

W_estW=estW(X_input,S_input,W_init);

figure;
subplot(121); imagesc(corr(X));
subplot(122); imagesc(corr(L));


figure;
subplot(121); imagesc(corr(methy));
subplot(122); imagesc(corr(expr));



figure;
subplot(121); imagesc(corr(X));
subplot(122);imagesc(corr(methy));



figure;
subplot(121);imagesc(corr(expr));
subplot(122);imagesc(corr(X*W_new));

runsim_estW_pca(702,5404,250,0.5,0,0,0,1)