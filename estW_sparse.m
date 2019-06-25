function W = estW_sparse(X,S,Wo,thr)
    K = 5;
    [n,m] = size(X);
    d = size(S,1);
%    W = randn(m,d);
%     W = Wo+randn(m,d);
%     W = double(Wo~=0);
%    W = double(Wo>0.25);
    W = Wo; 
    
    L=X*W;
    for k=1:K
        for i=1:d
            disp(i);
            options = optimoptions('fmincon','Display','iter');
            ind = W(:,i)~=0;
            fun=@(w)(-1*sum(([S(i,1:i-1) S(i,i+1:end)].*corr(X(:,ind)*w,[L(:,1:i-1) L(:,i+1:end)])))/d);
            W_raw = fmincon(fun,W(ind,i),[],[],[],[],[],[],[],options);
            W(ind,i) = wthresh(W_raw,'s',thr);
            L=X*W;
        end
    end
end