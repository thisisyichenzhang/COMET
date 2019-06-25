%% X: Data Matrix with number_of_rows=samples and number_of_cols=CpGs 
%% probe_gene_distance_matrix; rows= CpGs and cols =genes
%% W_initial = same dimension as probe_gene_distance_matrix
%% number of optimization iteration usually 4 is enough
%% gex = input network ; number of rows = number of cols = number of genes
function [W_learned_in_diff_iters, objective_function_in_diff_iters] = infer(X, probe_gene_distance_matrix,W_initial, number_of_optimization_iteration,gex)
        
    W_learned_in_diff_iters = cell(number_of_optimization_iteration, 1);

    objective_function_in_diff_iters = zeros(number_of_optimization_iteration+1, 1);
    
    
    W_learned = W_initial;
    number_of_genes = size(probe_gene_distance_matrix,2);

    L_learned=X*W_learned;

    objective_function_in_diff_iters(1)=mean(mean(diag(corr(corr(X*W_learned),gex))));

    nonlcon = @unitcircle;
    
    for opt_iter=1:1:number_of_optimization_iteration
        
        for i=1:1:number_of_genes
            disp(i);
            non_zeros=find(probe_gene_distance_matrix(:,i)~=0);
            options = optimoptions('fmincon','Display','iter');
            fun=@(w)(-1*sum(([gex(i,1:i-1) gex(i,i+1:end)].*corr(X(:,non_zeros)*w,[L_learned(:,1:i-1) L_learned(:,i+1:end)])))/number_of_genes);
            W_learned(non_zeros,i) = fmincon(fun,W_learned(non_zeros,i),[],[],[],[],[],[],nonlcon,options);
            L_learned=X*W_learned;
        end

        objective_function_in_diff_iters(opt_iter+1)=mean(mean(diag(corr(corr(X*W_learned),gex))));
        W_learned_in_diff_iters{opt_iter} = W_learned;

    end
end

