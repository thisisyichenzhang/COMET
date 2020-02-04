addpath('data/ROSMAP')
%addpath('../code')
%addpath(genpath('../../rosmapAD'))
%addpath(genpath('../../QUIC'))
%addpath('/nfs/37zfs1-packagces/Matlab/packages') % QUIUC and BIGQUIC path on orenthal

% parameters
chr = '1';
number_of_optimization_iterations=1;
upperlimit = 5000;
lowerlimit = 5000;
distance_based_initialization = 0;
%% Data Prepare
expressionAndPhenotype=load('expressionAndPhenotype');
methylationSNMnorm=load('methylationSNMnorm');
%methySNMtoGene100KB=load('methySNMtoGene100KB');
[~,meth_subjects,exp_subjects]=intersect(methylationSNMnorm.m.id,expressionAndPhenotype.data.projid);

probe_gene_distance_matrix = methySNMtoGene100KB.xToGene;
probes_on_this_chromosome = methylationSNMnorm.m.chr==str2double(chr);
probe_gene_distance_matrix=probe_gene_distance_matrix(probes_on_this_chromosome,:);

probe_gene_distance_matrix(probe_gene_distance_matrix<(-1*upperlimit) | probe_gene_distance_matrix>lowerlimit)=0;


X_transformed=zscore(methylationSNMnorm.m.data(:,meth_subjects)')';

X_transformed=X_transformed(probes_on_this_chromosome,:);
gene_symbols = expressionAndPhenotype.data.genes_symbols;

gene_symbols(sum(probe_gene_distance_matrix~=0,1)==0)=[];
probe_gene_distance_matrix(:,sum(probe_gene_distance_matrix~=0,1)==0)=[];
useless_CpGs=sum(probe_gene_distance_matrix~=0,2)==0;
probe_gene_distance_matrix(useless_CpGs,:)=[];
X_transformed(useless_CpGs,:)=[];
X_transformed=X_transformed';

gene_symbols(sum(probe_gene_distance_matrix~=0,1)==0)=[];

probe_gene_distance_matrix(:,sum(probe_gene_distance_matrix~=0,1)==0)=[];

% helper variables
number_of_genes = size(probe_gene_distance_matrix,2);
number_cpg_sites = size(probe_gene_distance_matrix,1);
number_of_subjects=size(X_transformed,1);

[a,b]=ismember(gene_symbols,expressionAndPhenotype.data.genes_symbols);
gex=expressionAndPhenotype.data.data(exp_subjects,b);

Gamma = distances_to_initial_weights_inverse_distance_proportional((probe_gene_distance_matrix~=0)',5000)';
if(distance_based_initialization==1)
   Gamma = distances_to_initial_weights_inverse_distance_proportional((probe_gene_distance_matrix)')';
end


W_initial = double(Gamma);


clear expressionAndPhenotype methylationSNMnorm methySNMtoGene100KB

file_name = strcat('ROSMAP_Final_upper_',num2str(upperlimit),'_lower_',num2str(lowerlimit),'_distance_based_',num2str(distance_based_initialization),'_chr1111_',num2str(chr));
    
[~,gex_netwok] = corr(gex);


% X_group=[];
% counter=1;
% W_initial=[];
% for gene_ind=1:1:number_of_genes
%    non_zeros=find(probe_gene_distance_matrix(:,gene_ind)~=0);
%     corr_local=corr(X_transformed(:,non_zeros));
%     [val,ind]=min(min(abs(corr_local)));
%     if(val<0.99)
%         [row,col]=ind2sub([size(corr_local,1) size(corr_local,2)],ind);
%         W_initial(counter,gene_ind)=1;
%         counter=counter+1;
%         W_initial(counter,gene_ind)=1;
%         counter=counter+1;
% 
%         X_group=[X_group X_transformed(:,non_zeros(row)) X_transformed(:,non_zeros(col))];
%     else
%         W_initial(counter,gene_ind)=1;
%         counter=counter+1;
%         X_group=[X_group X_transformed(:,gene_ind)];
%     end
% end
% X_transformed=X_group;

% X_transformed_group=[];
% virtual_cpg_counter=1;
% W_initial_grouped=[];
% for gene_ind=1:1:size(probe_gene_distance_matrix,2)
%     disp(gene_ind);
%     non_zeros=find(probe_gene_distance_matrix(:,gene_ind)~=0);
%     if(length(non_zeros)==1)
%        X_transformed_group=[X_transformed_group ]; 
%        W_initial_grouped(virtual_cpg_counter,gene_ind)=1;
%        virtual_cpg_counter = virtual_cpg_counter+1;
%        continue;
%     end
%     klist=1;%the number of clusters you want to try
%     myfunc = @(X,K)(kmeans(X, K,'Distance','correlation','Replicates',10));
%     
%     eva = evalclusters(X_transformed(:,non_zeros),myfunc,'CalinskiHarabasz','klist',klist);
%     optimal_k=1;
%     if(~isnan(eva.OptimalK))
%         optimal_k=eva.OptimalK;
%     end
%     classes=kmeans(X_transformed(:,non_zeros)',optimal_k,'Replicates',10);
%     for i=1:1:optimal_k
%        class_members = find(classes==i);
%        if(length(class_members)==1)
%           X_transformed_group=[X_transformed_group X_transformed(:,non_zeros(class_members))]; 
%        else
%           [pca_result,scores]=pca( X_transformed(:,non_zeros(class_members))','Centered',false);
%            X_transformed_group=[X_transformed_group pca_result(:,1)]; 
%        end
%        W_initial_grouped(virtual_cpg_counter,gene_ind)=1;
%        virtual_cpg_counter = virtual_cpg_counter+1;       
%     end
%     
% end

L_pca=[];
sign_val=[];
W_pca = zeros(size(Gamma,1),size(Gamma,2));
 for gene_ind=1:1:size(probe_gene_distance_matrix,2)
       disp(gene_ind)
      non_zeros=find(probe_gene_distance_matrix(:,gene_ind));
      if(length(non_zeros)==1)
          L_pca(:,gene_ind)=X_transformed(:,non_zeros);
        sign_val(gene_ind)=corr(L_pca(:,gene_ind),gex(:,gene_ind));
                  L_pca(:,gene_ind)=L_pca(:,gene_ind)*sign(corr(L_pca(:,gene_ind),gex(:,gene_ind)));

                  
          W_pca(non_zeros,gene_ind)=1;
          continue;
      end
      [pc_val,scores]=pca(X_transformed(:,non_zeros)','Centered',false);
        L_pca(:,gene_ind)=pc_val(:,1);
                sign_val(gene_ind)=corr(L_pca(:,gene_ind),gex(:,gene_ind));

        L_pca(:,gene_ind)=L_pca(:,gene_ind);%*sign(corr(L_pca(:,gene_ind),gex(:,gene_ind)));

 end
 
 
inds_perm = randperm(number_of_genes);
inds_perm = 1:1:number_of_genes;
beta=1;
[ W_in_diff_iters, objective_in_diff_iters] = infer(X_transformed,W_initial, W_initial, number_of_optimization_iterations,corr(gex(:,inds_perm)),beta,gex,L_pca);    

save(strcat(file_name,'.mat'),'-v7.3');



