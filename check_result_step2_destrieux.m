%% 

% This script will plot the simulated FC for HIV/HC.


clear all
close all
clc
figure_dir = '/gpfs/fs2/scratch/yzhuang4/CBIG-0.7.0-Wang2018_MFMem/result_figures';



%%  load the simulated FC for test HIV

% 
% h_output : simulated FC, only entries above main diagonal, in vector form
% CC_check: cross correlation of 2 FCs 
% CC_check = corr(atanh(FC_simR),atanh(y_FC));


% below is the simulation result calculated from the 
% best Para_E result calculated from 500 iterations, 9 random
% init ializeations.

 


load('agematched_1000simulation_testGrp_resLow_hiv_20191024.mat')
hiv.CC_check_test = CC_check;
hiv.h_output_test = h_output;

load('agematched_1000simulation_trainGrp_resLow_hiv_20191024.mat')
hiv.CC_check_train = CC_check;
hiv.h_output_train = h_output;


load('agematched_1000simulation_trainGrp_resLow_hiv_12wk_20191024.mat')
wk12.CC_check_train = CC_check;
wk12.h_output_train = h_output;


load('agematched_1000simulation_testGrp_resLow_hiv_12wk_20191024.mat')
wk12.CC_check_test = CC_check;
wk12.h_output_test = h_output;


load('agematched_1000simulation_trainGrp_resLow_hc_20191024.mat')
hc.CC_check_train = CC_check;
hc.h_output_train = h_output;


load('agematched_1000simulation_testGrp_resLow_hc_20191024.mat')
hc.CC_check_test = CC_check;
hc.h_output_test = h_output;


% 
% load('1000simulation_testGrp_resLow_hc_2019411.mat')
% CC_check_hc = CC_check;
% h_output_hc = h_output;

%%
numSimulation = 1000;

set(figure(1),'Position',[100 130 800 600])

subplot(312)
plot([1:numSimulation],hiv.CC_check_train,'bo-','markerfacecolor','b','markersize', 2)
% title('hiv training')
% xlabel('simulation number','FontSize',9)
% ylabel('Similarity','FontSize',9)
hold on
plot([1:numSimulation],hiv.CC_check_test,'ro-','markerfacecolor','r','markersize', 2)
ylim([0.32,0.62])
title(['hiv bsl mean similarity train : ', num2str(mean(hiv.CC_check_train),3), ' test : ', num2str(mean(hiv.CC_check_test),3)],'FontSize',15)
xlabel('simulation number','FontSize',15)
ylabel('Similarity','FontSize',15)
legend('train', 'test','Location','southeastoutside')


disp('hiv training goodness of fitting: ')
CC_avg = mean(hiv.CC_check_train)
disp('hiv test goodness of fitting: ')
CC_avg = mean(hiv.CC_check_test)


subplot(313)
plot([1:numSimulation],wk12.CC_check_train,'bo-','markerfacecolor','b','markersize', 2)
% title('hiv training')
% xlabel('simulation number','FontSize',9)
% ylabel('Similarity','FontSize',9)
hold on
plot([1:numSimulation],wk12.CC_check_test,'ro-','markerfacecolor','r','markersize', 2)
ylim([0.32,0.62])
title(['hiv wk12 mean similarity train : ', num2str(mean(wk12.CC_check_train),3), ' test : ', num2str(mean(wk12.CC_check_test),3)],'FontSize',15)
xlabel('simulation number','FontSize',15)
ylabel('Similarity','FontSize',15)
legend('train', 'test','Location','southeastoutside')


disp('hiv wk12 training goodness of fitting: ')
CC_avg = mean(wk12.CC_check_train)
disp('hiv wk12 test goodness of fitting: ')
CC_avg = mean(wk12.CC_check_test)

% plot result
%set(figure(2),'Position',[100 130 800 400])

subplot(311)
plot([1:numSimulation],hc.CC_check_train,'bo-','markerfacecolor','b','markersize', 2)
% title('hc training')
% xlabel('simulation number','FontSize',9)
% ylabel('Similarity','FontSize',9)

%subplot(212)
hold on
ylim([0.32,0.62])
plot([1:numSimulation],hc.CC_check_test,'ro-','markerfacecolor','r','markersize', 2)
title(['hc mean similarity train : ', num2str(mean(hc.CC_check_train),3), ' test : ', num2str(mean(hc.CC_check_test),3)],'FontSize',15)
xlabel('simulation number','FontSize',15)
ylabel('Similarity','FontSize',15)
legend('train', 'test','Location','southeastoutside')


disp('hc training goodness of fitting: ')
CC_avg = mean(hc.CC_check_train)
disp('hc test goodness of fitting: ')
CC_avg = mean(hc.CC_check_test)


%saveas(figure(1), [figure_dir, '/agematched_step2_1000sim_result_rmltnp_destrieux.fig'])
%saveas(figure(1), [figure_dir, '/agematched_step2_1000sim_result_rmltnp_destrieux.pdf'])


%% calculate the mean zFC from 1000 simulation

%fisher-z transform the FC for HIV, HC
for i = 1: numSimulation
    hz_hc_test(i,:) = atanh(hc.h_output_test(i,:));
    hz_hc_train(i,:) = atanh(hc.h_output_train(i,:));
    hz_hiv_test(i,:) = atanh(hiv.h_output_test(i,:));
    hz_hiv_train(i,:) = atanh(hiv.h_output_train(i,:));
    hz_wk12_test(i,:) = atanh(wk12.h_output_test(i,:));
    hz_wk12_train(i,:) = atanh(wk12.h_output_train(i,:));
end

hc.h_avg = mean(hz_hc_test,1);
set(figure(2),'Position',[100 130 800 600])

clims = [0,1];

hc.FCsim = triuvec2matrix(148, hc.h_avg);
figure(2),subplot(332),
imagesc(hc.FCsim,clims)
title('hc simulated fc')
colorbar()

hiv.h_avg = mean(hz_hiv_test,1);

hiv.FCsim = triuvec2matrix(148, hiv.h_avg);
figure(2),subplot(335),
imagesc(hiv.FCsim,clims)
title('hiv bsl simulated fc')
colorbar()

wk12.h_avg = mean(hz_wk12_test,1);

wk12.FCsim = triuvec2matrix(148, wk12.h_avg);
figure(2),subplot(338),
imagesc(wk12.FCsim,clims)
title('hiv wk12 simulated fc')
colorbar()


% load empirical test FC
load('train_test_dataset_agematch_destrieux_1004.mat')

hc_fc_test = hc_fc_test - diag(diag(hc_fc_test));
figure(2),subplot(331),
imagesc(hc_fc_test,clims), colorbar
title('hc empirical fc')

hiv_fc_test = hiv_fc_test - diag(diag(hiv_fc_test));
figure(2), subplot(334)
imagesc(hiv_fc_test, clims), colorbar
title('hiv bsl empirical fc')


wk12_fc_test = hiv_fc_test_12wk - diag(diag(hiv_fc_test_12wk));
figure(2), subplot(337)
imagesc(hiv_fc_test, clims), colorbar
title('hiv wk12 empirical fc')


%
num_rois = 148;
tri_u_mask = ones(num_rois);
tri_u_mask = triu(tri_u_mask,1);
%figure,imagesc(tri_u_mask)
%figure,imagesc(hc_fc_test)
a = hc_fc_test(tri_u_mask == 1);

% % add linear fit line to the plot
% beta_hc = h_avg_hc'\a;
% y2_hc = a*beta_hc;

figure(2),
subplot(333)
scatter(a, hc.h_avg,5, 'filled');
hold on
ylim([-0.2,1])
xlim([-0.4,1])
disp('HC similarity : ')
[r,p] = corr(a,hc.h_avg')
% plot(a,y2_hc,'r-')
l1 = lsline;
set(l1, 'lineWidth', 2, 'Color', 'r')
hold off
xlabel('zFC empirical')
ylabel('zFC simulated')
title(['Correlation r=', num2str(corr(a,hc.h_avg'),3)])

b = hiv_fc_test(tri_u_mask == 1);
c = hiv_fc_test_12wk(tri_u_mask == 1);
% % add linear fit line to the plot
% beta_hiv = h_avg_hiv'\b;
% y2_hiv = b*beta_hiv;

figure(2),
subplot(336)
scatter(b, hiv.h_avg,5, 'filled');
hold on
ylim([-0.2,1])
xlim([-0.4,1])
disp('HIV bsl similarity : ')

[r,p] = corr(b,hiv.h_avg')
% plot(b,y2_hiv,'r-');
l2 = lsline;
set(l2, 'lineWidth', 2, 'Color', 'r')
hold off
xlabel('zFC empirical')
ylabel('zFC simulated')
title(['Correlation r=', num2str(corr(b,hiv.h_avg'),3)])

figure(2),
subplot(339)
scatter(c, wk12.h_avg, 5, 'filled');
hold on
ylim([-0.2,1])
xlim([-0.4,1])
disp('HIV wk12 similarity : ')

[r,p] = corr(c,wk12.h_avg')
l2 = lsline;
set(l2, 'lineWidth', 2, 'Color', 'r')
hold off
xlabel('zFC empirical')
ylabel('zFC simulated')
title(['Correlation r=', num2str(corr(c,wk12.h_avg'),3)])

%saveas(figure(2), [figure_dir, '/agematched_step2_avg_FC_sim_rmltnp_destrieux.fig'])
%saveas(figure(2), [figure_dir, '/agematched_step2_avg_FC_sim_rmltnp_destrieux.pdf'])


%% check the meanFC distribution of HIV and HC



%% use simulated FC, randomly choose 30 hiv and 30 hc from test dataset.

rng(1)
num_subj = 50;
num_rois = 148;

subjidx_hiv = randi(numSimulation, 1, num_subj);
subjidx_hc = randi(numSimulation, 1, num_subj);
subjidx_wk12 = randi(numSimulation, 1, num_subj);

fcm_sim_hiv = zeros(num_rois,num_rois,num_subj);
fcm_sim_hc = zeros(num_rois,num_rois,num_subj);
fcm_sim_wk12 = zeros(num_rois,num_rois,num_subj);

for idx = 1:num_subj
    fcm_sim_hiv(:,:, idx) = triuvec2matrix(68, hiv.h_output_test(idx,:)) ;
    fcm_sim_hc(:,:, idx) = triuvec2matrix(68, hc.h_output_test(idx,:)) ;
    fcm_sim_wk12(:,:, idx) = triuvec2matrix(68, wk12.h_output_test(idx,:)) ;

end




%% compare edges for simulated FC
clear fctemp
clear s
clear p

p = zeros(num_rois);
grp = [repmat({'hc'},1,num_subj),repmat({'hiv'},1,num_subj),repmat({'wk12'},1,num_subj)];

% for idx_i = 1:10
%     for idx_j = idx_i:10        
%         fctemp = [squeeze(fcm_sim_hc(idx_i,idx_j,:)),squeeze(fcm_sim_hiv(idx_i,idx_j,:)),squeeze(fcm_sim_wk12(idx_i,idx_j,:))];
%         [p(idx_i,idx_j), s{idx_i, idx_j}] = anova1(fctemp, grp,'on');
%     end
% end

% use p = 0.05 as threshold
p_thr = 0.001;

% FDR correction
[h_fdr,p_fdr] = fdr_bh(p, p_thr);

figure(3),
imagesc(h_fdr),colorbar
title('two sample t-test HIV vs. HC simulated FC')

%% calculate the empirical sc-fc similarity in test dataset
disp('HC sc-fc empirical similarity: ')
num_rois = 148;
tri_u_mask = ones(num_rois);
tri_u_mask = triu(tri_u_mask,1);
[r,p] = corr(hc_fc_test(tri_u_mask == 1),hc_sc_test(tri_u_mask == 1) )


disp('HIV BSL sc-fc empirical similarity: ')
[r,p] = corr(hiv_fc_test(tri_u_mask == 1),hiv_sc_test(tri_u_mask == 1) )


disp('HIV wk12 sc-fc empirical similarity: ')
[r,p] = corr(hiv_fc_test_12wk(tri_u_mask == 1),hiv_sc_test_12wk(tri_u_mask == 1) )




