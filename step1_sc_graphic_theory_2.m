clear all
close all
clc

%% 1) load age matched data

p = genpath('/Users/yuchuan/Box Sync/matlab_toolbox');
addpath(p);
p = genpath('/Users/yuchuan/Box Sync/Connectome/scripts');
addpath(p);

cd('/Users/yuchuan/Box Sync/Connectome/')
load age_matched_data.mat

age = [subj_list_agematch{:,2}];

%% 2) calculate nodal graphic theroy parameters, correlate with w, I

visit_str = {'BSL', 'Wk12'};
num_visit = length(visit_str);


% number of random network generate
max_iter = 5;

% define sparsity for SC
sparsity = linspace(0.1,0.4,10);

num_spar = length(sparsity);


% ignore FC with 1, and 10
% 1 - lh_lateral ventricle
% 10 - rh_lateral_ventricle
%roi_index = [2:9,11:87];
roi_index = [20:87];


% clustering coefficient
% degree
% local efficiency
clus_coef = zeros([num_subject,num_visit,num_spar,length(roi_index)]);
degree = zeros([num_subject,num_visit,num_spar,length(roi_index)]);
local_efficiency = zeros([num_subject,num_visit,num_spar,length(roi_index)]);
sc_feature = 2; % sc feature 1: max fa
                % sc feature 2: mean fa; 
                % sc feature 3: max md; 
                % sc feature 4: mean md; 
                % sc feature 5: count   
                % sc feature 6: fiber length 
                % sc feature 7: connective volumes 
                % sc feature 8: connective volume ratio

                
% for each entry of the matrix, if there
% were less than 70% subjects with fibers in the
% entry, it was set to zero.
for i = 1:87
    for j = i:87
        %if  nnz(squeeze(scm(i,j,:,1,sc_feature))) < 46
        if  nnz(squeeze(scm(i,j,:,1,sc_feature))) < round(92*0.7)
            [i,j]
            scm(i,j,:,:,:) = 0;
        end
    end
end


for idx_subj = 1:num_subject
    idx_subj
    for idx_visit = 1:num_visit
% for idx_subj = 1:4
%     for idx_visit = 1:1
        
        %temp_b = abs(temp);
        temp=scm(roi_index,roi_index,idx_subj,idx_visit,sc_feature);
        temp_b = temp + triu(temp, -1).';
        if max(max(temp_b)) ~= 0;
            temp_b = weight_conversion(temp_b,'normalize');
            for idx_sparsity = 1:length(sparsity)

                temp_bin = threshold_proportional(temp_b, sparsity(idx_sparsity));
                clus_coef(idx_subj,idx_visit,idx_sparsity,:) = clustering_coef_wu(temp_bin);
                local_efficiency(idx_subj,idx_visit,idx_sparsity,:) = efficiency_wei(temp_bin,2);
                degree(idx_subj,idx_visit,idx_sparsity,:) = degrees_und(temp_bin);
        
%                 SPL = distance_wei_floyd(temp_bin,'inv');
%                 [lambda(idx_subj,idx_visit,idx_sparsity), efficiency(idx_subj,idx_visit,idx_sparsity) ] = charpath(SPL);

            end
        end
    end
end

%% 

node_idx = 7;
clus_coef_nodal = clus_coef(:,:,:,node_idx);
clus_coef_nodal(clus_coef_nodal == 0) = nan;


degree_nodal = degree(:,:,:,node_idx);
degree_nodal(degree_nodal == 0) = nan;

efficiency_nodal = local_efficiency(:,:,:,node_idx);

efficiency_nodal(efficiency_nodal == 0) = nan;
% betweenness_net(betweenness_net == 0) = nan;
% eigenvector_centrality_net(eigenvector_centrality_net == 0) = nan;
% assortativity(assortativity == 0) = nan;

n_sub = num_subject;
n_hiv = num_hiv;
n_hc = n_sub - n_hiv;

% smallworldness = sigma;


clus_coef_1d_3groups = [reshape(clus_coef_nodal(1:n_hiv,1,:),1,[]),reshape(clus_coef_nodal(1:n_hiv,2,:),1,[]),reshape(clus_coef_nodal(n_hiv+1 : end ,1,:),1,[]),];
degrees_1d_3groups = [reshape(degree_nodal(1:n_hiv,1,:),1,[]),reshape(degree_nodal(1:n_hiv,2,:),1,[]),reshape(degree_nodal(n_hiv+1 : end ,1,:),1,[]),];
efficiency_1d_3groups = [reshape(efficiency_nodal(1:n_hiv,1,:),1,[]),reshape(efficiency_nodal(1:n_hiv,2,:),1,[]),reshape(efficiency_nodal(n_hiv+1 : end ,1,:),1,[]),];
% clus_coef_norm_1d_3groups = [reshape(clus_coef_norm(1:n_hiv,1,:),1,[]),reshape(clus_coef_norm(1:n_hiv,2,:),1,[]),reshape(clus_coef_norm(n_hiv+1 : end ,1,:),1,[]),];
% lambda_norm_1d_3groups = [reshape(lambda_norm(1:n_hiv,1,:),1,[]),reshape(lambda_norm(1:n_hiv,2,:),1,[]),reshape(lambda_norm(n_hiv+1 : end ,1,:),1,[]),];
% sigma_1d_3groups = [reshape(sigma(1:n_hiv,1,:),1,[]),reshape(sigma(1:n_hiv,2,:),1,[]),reshape(sigma(n_hiv+1 : end ,1,:),1,[]),];
% modularity_1d_3groups = [reshape(modularity(1:n_hiv,1,:),1,[]),reshape(modularity(1:n_hiv,2,:),1,[]),reshape(modularity(n_hiv+1 : end ,1,:),1,[]),];
% eigenvector_centrality_1d_3groups = [reshape(eigenvector_centrality_net(1:n_hiv,1,:),1,[]),reshape(eigenvector_centrality_net(1:n_hiv,2,:),1,[]),reshape(eigenvector_centrality_net(n_hiv+1 : end ,1,:),1,[]),];
% betweenness_1d_3groups = [reshape(betweenness_net(1:n_hiv,1,:),1,[]),reshape(betweenness_net(1:n_hiv,2,:),1,[]),reshape(betweenness_net(n_hiv+1 : end ,1,:),1,[]),];
% assortativity_1d_3groups = [reshape(assortativity(1:n_hiv,1,:),1,[]),reshape(assortativity(1:n_hiv,2,:),1,[]),reshape(assortativity(n_hiv+1 : end ,1,:),1,[]),];
%group_vec = [repmat({'hiv bsl'}, 1, n_hiv), repmat({'hiv 12wk'}, 1, n_hiv), repmat({'hc'}, 1, n_hc)];
group_vec = [repmat({'hiv bsl'}, 1, n_hiv * num_spar), repmat({'hiv 12wk'}, 1, n_hiv * num_spar), repmat({'hc'}, 1, n_hc * num_spar)];
group_1d = group_vec;

sparsity_1d = [reshape( repmat(sparsity, n_hiv ,1 ), 1, [] ),reshape( repmat(sparsity, n_hiv ,1 ), 1, [] ),reshape( repmat(sparsity, n_hc ,1 ), 1, [] )];

close all

figure(1)
g1 = gramm('x', sparsity_1d, 'y', clus_coef_1d_3groups, 'color', group_1d);
%g1.geom_point('alpha',0.1);
g1.stat_boxplot('notch', true);

g1.set_title('Cluster Coefficient');
g1.set_names('x','Sparsity','y','Cluster Coefficient');
g1.set_text_options('label_scaling',2,'legend_scaling',2);
g1.axe_property('FontSize', 20);
g1.draw();

figure(2)
g2 = gramm('x', sparsity_1d, 'y', degrees_1d_3groups, 'color', group_1d);
g2.stat_boxplot('notch', true);
g2.set_title('degree');
g2.set_names('x','Sparsity','y','degree');
g2.set_text_options('label_scaling',2,'legend_scaling',2);
g2.axe_property('FontSize', 20);
g2.draw();

figure(3)
g3 = gramm('x', sparsity_1d, 'y', efficiency_1d_3groups, 'color', group_1d);
g3.stat_boxplot();
g3.set_title('Efficiency - SC');
g3.set_names('x','Sparsity','y','Efficiency');
g3.set_text_options('label_scaling',2,'legend_scaling',2);
g3.axe_property('FontSize', 20);
g3.draw();



%%
hiv_clus_coef = squeeze(mean(clus_coef(1:num_hiv,1,6,:),1));
hiv_degree= squeeze(mean(degree(1:num_hiv,1,6,:),1));
hiv_efficiency_local = squeeze(mean(local_efficiency(1:num_hiv,1,6,:),1));

hc_clus_coef = squeeze(mean(clus_coef(47:end,1,6,:),1));
hc_degree = squeeze(mean(degree(47:end,1,6,:),1));
hc_efficiency_local = squeeze(mean(local_efficiency(47:end,1,6,:),1));


wk12_clus_coef = squeeze(mean(clus_coef(1:num_hiv,2,6,:),1));
wk12_degree = squeeze(mean(degree(1:num_hiv,2,6,:),1));
wk12_efficiency_local = squeeze(mean(local_efficiency(1:num_hiv,2,6,:),1));



%% 2) calculate graphic theory metrics
% 2.2) binary% clustering coefficient
% degree
% local efficiency
clus_coef = zeros([num_subject,num_visit,num_spar,length(roi_index)]);
degree = zeros([num_subject,num_visit,num_spar,length(roi_index)]);
local_efficiency = zeros([num_subject,num_visit,num_spar,length(roi_index)]);


for idx_subj = 1:num_subject
    idx_subj
    for idx_visit = 1:num_visit
% for idx_subj = 1:4
%     for idx_visit = 1:1
        temp=scm(roi_index,roi_index,idx_subj,idx_visit,sc_feature);
        temp_b = temp + triu(temp, -1).';

        if max(max(temp_b)) ~= 0;
            temp_b = weight_conversion(temp_b,'normalize');
            
            for idx_sparsity = 1:length(sparsity)
                
                temp_bin = threshold_proportional(temp_b, sparsity(idx_sparsity));
                temp_bin = weight_conversion(temp_bin,'binarize');

                clus_coef(idx_subj,idx_visit,idx_sparsity,:) = clustering_coef_bu(temp_bin);
                degree(idx_subj,idx_visit,idx_sparsity,:) = degrees_und(temp_bin);
                local_efficiency(idx_subj,idx_visit,idx_sparsity,:)  = efficiency_bin(temp_bin, 1);
%                 SPL = distance_bin(temp_bin);
%                 [lambda(idx_subj,idx_visit,idx_sparsity), efficiency(idx_subj,idx_visit,idx_sparsity) ] = charpath(SPL);
%         
%                 if max(max(temp_bin ~= 0))
%                     [Ci(idx_subj,idx_visit,idx_sparsity,:), Q(idx_subj,idx_visit,idx_sparsity,:) ] = modularity_und(temp_bin);
%                 end
            end
        end
    end
end

%% 







