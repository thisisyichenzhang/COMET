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
    X(:,ind) = sig*ones(1,b) + 0.1*randn(n,b);
    W(ind,i) = rand(b,1);
end

L = X*W;
S = corr(L);
Wl = estW(X,S,W);

function W = estW(X,S,Wo)
    K = 5;
    [n,m] = size(X);
    d = size(S,1);
%    W = randn(m,d);
%     W = Wo+randn(m,d);
     W = double(Wo~=0);
    W = double(Wo>0.25);
    L=X*W;
    for k=1:K
        for i=1:d
            disp(i);
            options = optimoptions('fmincon','Display','iter');
            ind = W(:,i)~=0;
            fun=@(w)(-1*sum(([S(i,1:i-1) S(i,i+1:end)].*corr(X(:,ind)*w,[L(:,1:i-1) L(:,i+1:end)])))/d);
            W(ind,i) = fmincon(fun,W(ind,i),[],[],[],[],[],[],[],options);
            L=X*W;
        end
    end
end

