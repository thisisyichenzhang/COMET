%indNZ = find(W(:)~=0);
%Ws = zeros(length(indNZ),1);
%Wls = zeros(length(indNZ),1);
%group_corr = zeros(d,1);
%for i = 1:d
%    ind = (1+(i-1)*b):(i*b);
%    Ws(ind) = zscore(W(indNZ(ind)));
%    Wls(ind) = zscore(Wl(indNZ(ind)));
%    group_corr(i) = corr(Ws(ind), Wls(ind));
%end

%figure; 
%subplot(1,2,1);
%plot(zscore(W(W~=0)));
%hold on; plot(zscore(Wl(W~=0)));
%title("Ture zscore of W and estimated W");
%subplot(1,2,2);
%plot(Ws);
%hold on; plot(Wls);
%title("Ture zscore of Group-wise W and estimated W");


%figure; 
%subplot(1,2,1); stem(diag(corr(W,Wl))); title("Correlations W & W0 ");
%subplot(1,2,2); stem(group_corr); title(" Group-wise correlations W & W0 ");
%save('dfdf.mat')


% Visualisation
W = W_cell{1};
X = X_cell{1};
L = X*W;
S = corr(L);
We = We_cell{1};
% Estimated L, and S
Le = X * We;
Se = corr(Le);
    
    
figure;
subplot(1,2,1);
imagesc(S.*~eye(d));
title("gene-wise corr (L)");
subplot(1,2,2);
imagesc(Se.*~eye(d));
title("Gene wise corr (L_{e})");


figure; plot(zscore(W(W~=0)));
hold on; plot(zscore(We(W~=0)));

figure;
subplot(1,2,1);
imagesc(L);
title("L");
subplot(1,2,2);
imagesc(Le);
title("L_{e}");