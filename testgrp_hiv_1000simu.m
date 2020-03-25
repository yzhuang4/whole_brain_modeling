
clear all
close all
clc

%% setup directory for lib, data, save
main_dir = pwd;
data_dir = fullfile(main_dir,'data'); 
save_dir = pwd;

cd('..')
high_dir = pwd;
lib_dir = fullfile(high_dir, 'lib');

addpath(lib_dir);



%% load model parameter
%parameter_file_name = 'Estimated_Parameter_hiv_2019411.mat';
%parameter_file_name = '500iters_final_estimate_Paras_agematched_0708.mat' ;
parameter_file_name = '500iters_final_estimate_Paras_agematched_1003.mat' ;

load ([data_dir '/' parameter_file_name],'hiv');

Para_E = hiv.Para_E_avg;



%% load FC, SC
% FCSC_file_name = 'FCSC_Desikan68_Raphael_Wang.mat';
% load([data_dir '/' FCSC_file_name]);
% 
% SC = SC_test;
% FC = FC_test;
% 
% 
% %scaling the SC
% SC = SC./max(max(SC)).*0.2;
% 
% %prepare FC: use the entries above main diagonal
% FC_mask = tril(ones(size(FC,1),size(FC,1)),0);
% y = FC(~FC_mask);



%load(fullfile(data_dir,'train_test_dataset_agematch_0625.mat'))
load(fullfile(data_dir,'train_test_dataset_agematch_0927.mat'))

FC = hiv_fc_test;
FC = FC ./ max(max(FC));
SC = hiv_sc_test; 
SC = SC./max(max(SC)).*0.2;
FC_mask = tril(ones(size(FC,1),size(FC,1)),0);
y = FC(~FC_mask); %use the elements above the maiin diagnal, y becomes a vector {samples x 1} 


%% begin simulation


%funcP = @(Para_E,Nstate, BOLD_d) CBIG_MFMem_rfMRI_nsolver_eul_sto_resLH(Para_E,SC,y,FC_mask,Nstate,14.4,0.72,0);
funcP = @(Para_E,Nstate, BOLD_d) CBIG_MFMem_rfMRI_nsolver_eul_sto_resLH(Para_E,SC,y,FC_mask,Nstate,7,2,0);
  
numSimulation = 1000;

for i = 1:numSimulation

disp(['Sim:' num2str(i)]);

Nstate = rng;
[h_output(i,:), CC_check(i), BOLD_d(i,:,:)] = funcP(Para_E,Nstate);

CC_check(i)

end

%% plot result
set(figure,'Position',[100 130 600 400],'Color','w')
plot([1:length(CC_check)],CC_check,'ko-','markerfacecolor','k')

%plot([1:numSimulation],CC_check,'ko-','markerfacecolor','k')
xlabel('simulation number','FontSize',9)
ylabel('Similarity','FontSize',9)

%% save result
saved_date = fix(clock);

save( [save_dir '/agematched_1000simulation_testGrp_resLow_hiv_' num2str(saved_date(1)) num2str(saved_date(2)) num2str(saved_date(3))],'CC_check','h_output','BOLD_d');
rmpath(lib_dir);
