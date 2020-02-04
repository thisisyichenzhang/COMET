%% Run Cross-validation on real data 
% @X_input: samples by features 
% @S_input: co-expression network
% @W_init: initilization of W
% @nfolfs: number of folds for cross-validatin 
% @seed: random seed for reproducibility

% Author: Yichen Zhang (yichen.zhang@stat.ubc.ca)
% Last Update: Aug 13th 2019.
function cv_estW_smooth(X_input,S_input,W_init,nfolds,seed)
        % Normalization of X
        X_input = normalize(X_input);
        n = size(X_input,1);
        % The range of parameters for tuning
        alphas = [0:0.01:0.1,0.2,0.5,0.8,1];
        % Set the seed for reproducibility.
        rng(seed);
        
        % Perform Cross-validation
        mean_criterion_byalpha = zeros(length(alphas),1);
        % Cross-validation Partition
        CVO = cvpartition(n,'k',nfolds);
        criterion_byfold = zeros(CVO.NumTestSets,1);

        for k = 1:length(alphas)
            for i = 1:CVO.NumTestSets
                trIdx = CVO.training(i);
                teIdx = CVO.test(i);
                [~,criterion] = estW_smooth(X_input(trIdx,:),X_input(teIdx,:),S_input,W_init,alphas(k));
                criterion_byfold(i) = criterion;  
            end
            mean_criterion_byalpha(k) = mean(criterion_byfold);

        end
    [~, Idx_max_criteirion] = max(mean_criterion_byalpha);
    best_alpha = alphas(Idx_max_criteirion);

    fprintf('Best alpha is found to be %0.2f\n', best_alpha);  

    % Runing on the full data with best tuning parameter found by CV
    [Wl_cv,~] = estW_smooth(X_input,X_input,S_input,W_init,best_alpha);
    
    % Saving results
    fname = sprintf('myfilename');
    save(fname,'Wl_cv','best_alpha','-v7.3');
    disp('Results saved.');
  
end
