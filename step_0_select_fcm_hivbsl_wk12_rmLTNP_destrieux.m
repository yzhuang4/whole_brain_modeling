%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step0_1_select_scmfcm.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Age matched group, remove some severe headmotion subjects, rerun after
%slice timing correction
% hiv bsl, hiv wk12

clear all
close all
clc
%% load sc and fc
cd('/Users/yuchuan/Box Sync/Connectome/')
load scm_destrieux.mat
load fcm_destrieux.mat
load demographic.mat
age = [demo_subjID{:,2}];

load('train_test_dataset_agematch_0927_all.mat', 'train_hc','train_hiv','train_wk12','test_hc','test_hiv','test_wk12');

%% get the exclusion list, 1 - keep, 0 - exclude
filename = '/Users/yuchuan/Box Sync/Connectome/scripts/agematched_rerun/age_match_exclude_LTNP.txt';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, '%d', 'ReturnOnError', false);
fclose(fileID);
exclusion_idx = [dataArray{:}];
%% 
count = 1;
for i = 1:length(exclusion_idx)
    if exclusion_idx(i) == 1
        subj_list_agematch{count,1} = demo_subjID{i,1};
        subj_list_agematch{count,2} = demo_subjID{i,2};
        subj_list_agematch{count,3} = demo_subjID{i,3};
        scm(:,:,count,:,:) = connectome_destrieux(:,:,i,:,:);
        fcm(:,:,count,:) = fcm_destrieux_2(:,:,i,:);
        zfcm(:,:,count,:) = fcm_z_destrieux_2(:,:,i,:);
        count = count + 1;
    end
end
%% check the subj_list_agematch cell array, change the number of subjects in each group
num_hiv = 42;
num_hc = 46;
num_subject = 88;
%% 
%% Split into training and testing dataset.
% ramdomly split HIV+ and HC subjects into half, test the age difference 

age = [subj_list_agematch{:,2}];
idx_hiv = 1:num_hiv;
idx_wk12 = idx_hiv;
idx_hc = num_hiv+1 : num_subject;
% if the subject's dti/fmri data is missing, make the index = nan

for i = 1:num_hiv
    if (max(max(scm(:,:,i,1,2)))==0) || (max(max(fcm(:,:,i,1))) == 0)
        idx_hiv(i) = nan
    end
end

for i = 1:num_hiv
    if (max(max(scm(:,:,i,2,2)))==0) || (max(max(fcm(:,:,i,2))) == 0)
        idx_wk12(i) = nan
    end
end

for i = 1 : num_hc
    if (max(max(scm(:,:,num_hiv+i,1,2)))==0) || (max(max(fcm(:,:,num_hiv+i,1))) == 0)
        idx_hc(i) = nan
    end
end

idx_hiv_rm = idx_hiv(~isnan(idx_hiv));
idx_wk12_rm = idx_hiv(~isnan(idx_wk12));
idx_hc_rm = idx_hc(~isnan(idx_hc));


num_hiv = length(idx_hiv_rm)
num_wk12 = length(idx_wk12_rm)
num_hc = length(idx_hc_rm)

num_trainhiv = round(num_hiv / 2);
num_testhiv = num_hiv - num_trainhiv;

num_trainwk12 = round(num_wk12 / 2);
num_testwk12 = num_wk12 - num_trainwk12;

num_trainhc = round(num_hc / 2);
num_testhc = num_hc - num_trainhc;

%%
% h1 = 1; 
% h2 = 1; 
% h3 = 1; 
% 
% %counter = 1;
% 
% %while h1 + h2 >= 1
% 
% for iter = 1:5000;
%     
%     randidx_hiv = randperm(numel(idx_hiv_rm));
%     tmp_train_hiv(iter,:) = idx_hiv_rm(randidx_hiv(1:num_trainhiv));
%     tmp_test_hiv(iter,:) = idx_hiv_rm(randidx_hiv(num_trainhiv + 1 : end));
% 
%     randidx_wk12 = randperm(numel(idx_wk12_rm));
%     tmp_train_wk12(iter,:) = idx_wk12_rm(randidx_wk12(1:num_trainwk12));
%     tmp_test_wk12(iter,:) = idx_wk12_rm(randidx_wk12(num_trainwk12 + 1 : end));
% 
%     randidx_hc = randperm(numel(idx_hc_rm));
%     tmp_train_hc(iter,:) = idx_hc_rm(randidx_hc(1:num_trainhc));
%     tmp_test_hc(iter,:) = idx_hc_rm(randidx_hc(num_trainhc+1 : end));
% 
%     [h1,p1] = ttest2(age(tmp_train_hiv(iter,:)), age(tmp_test_hiv(iter,:)));
%     [h3,p3] = ttest2(age(tmp_train_hc(iter,:)), age(tmp_test_hc(iter,:)));
%     [h2,p2] = ttest2(age(tmp_train_wk12(iter,:)), age(tmp_test_wk12(iter,:)));
%     
%     [h4,p4] = ttest2(age(tmp_train_hiv(iter,:)), age(tmp_train_hc(iter,:)));
%     [h5,p5] = ttest2(age(tmp_train_wk12(iter,:)), age(tmp_train_hc(iter,:)));
%     [h6,p6] = ttest2(age(tmp_train_hiv(iter,:)), age(tmp_train_wk12(iter,:)));
%     
%     [h7,p7] = ttest2(age(tmp_test_hiv(iter,:)), age(tmp_test_hc(iter,:)));
%     [h8,p8] = ttest2(age(tmp_test_wk12(iter,:)), age(tmp_test_hc(iter,:)));
%     [h9,p9] = ttest2(age(tmp_test_hiv(iter,:)), age(tmp_test_wk12(iter,:)));
%    
%     %counter = counter + 1;
%     h(iter) =h1+h2+h3+h4+h5+h6+h7+h8+h9;
%     
%     diff_age(iter,1) = mean(age(tmp_train_hiv(iter,:))) - mean(age(tmp_test_hiv(iter,:)));
%     diff_age(iter,3) = mean(age(tmp_train_hc(iter,:))) - mean(age(tmp_test_hc(iter,:)));
%     diff_age(iter,2) = mean(age(tmp_train_wk12(iter,:))) - mean(age(tmp_test_wk12(iter,:)));
%     diff_age(iter,4) = abs(diff_age(iter,1)) + abs(diff_age(iter,2))+ abs(diff_age(iter,3));
%     
% end
% 
% %%
% [age_min, idx] = min(diff_age(:,4))
% train_hiv = tmp_train_hiv(idx,:);
% test_hiv = tmp_test_hiv(idx,:);
% train_wk12 = tmp_train_wk12(idx,:);
% test_wk12 = tmp_test_wk12(idx,:);
% train_hc = tmp_train_hc(idx,:);
% test_hc = tmp_test_hc(idx,:);
% 
% [h1,p1] = ttest2(age(train_hiv), age(test_hiv))
% [h2,p2] = ttest2(age(train_wk12), age(test_wk12))
% [h3,p3] = ttest2(age(train_hc), age(test_hc))
% 
% std(age(train_hiv))
% std(age(test_hiv))
% std(age(train_wk12))
% std(age(test_wk12))
% std(age(train_hc))
% std(age(test_hc))

%% Use the index to calculate the mean SC_train, mean FC_train, mean SC_test, FC_test
% only extract the SC including the cortical regions 

num_rois = 148;

hiv_sc_train = zeros(num_rois);
hiv_fc_train = zeros(num_rois);
hiv_sc_test = zeros(num_rois);
hiv_fc_test = zeros(num_rois);

wk12_sc_train = zeros(num_rois);
wk12_fc_train = zeros(num_rois);
wk12_sc_test = zeros(num_rois);
wk12_fc_test = zeros(num_rois);

hc_sc_train = zeros(num_rois);
hc_fc_train = zeros(num_rois);
hc_sc_test = zeros(num_rois);
hc_fc_test = zeros(num_rois);



for i = 1:num_rois
    for j = 1:num_rois
        
        % To create the group-averaged SC matrix, for each entry of the matrix, if there
        % were less than 50% of the training (or test) subjects with fibers in the
        % entry, it was set to zero.
        
        % due to the low number of subjects in each cohort/train-test dataset 
        % combine all cohort for sc matrix check 
        
        if  nnz(squeeze(scm(19+i,19+j,:,1,5))) >= 46
            hiv_sc_train(i,j) = mean(scm(19+i,19+j,train_hiv,1,5));
            hiv_sc_test(i,j) = mean(scm(19+i,19+j,test_hiv,1,5));
            hc_sc_train(i,j) = mean(scm(19+i,19+j,train_hc,1,5));
            hc_sc_test(i,j) = mean(scm(19+i,19+j,test_hc,1,5));
            hiv_sc_train_12wk(i,j) = mean(scm(19+i,19+j,train_wk12,2,5));
            hiv_sc_test_12wk(i,j) = mean(scm(19+i,19+j,test_wk12,2,5));       
        else        
            hiv_sc_train(i,j) = 0;
            hiv_sc_test(i,j) = 0;
            hc_sc_train(i,j) = 0;
            hc_sc_test(i,j) = 0;
            hiv_sc_train_12wk(i,j) = 0;
            hiv_sc_test_12wk(i,j) = 0;
        end
        
        hiv_fc_train(i,j) = mean(fcm(19+i,19+j,train_hiv,1));
        hiv_fc_test(i,j) = mean(fcm(19+i,19+j,test_hiv,1));
        hc_fc_train(i,j) = mean(fcm(19+i,19+j,train_hc,1));
        hc_fc_test(i,j) = mean(fcm(19+i,19+j,test_hc,1));
        hiv_fc_train_12wk(i,j) = mean(fcm(19+i,19+j,train_wk12,2));
        hiv_fc_test_12wk(i,j) = mean(fcm(19+i,19+j,test_wk12,2));

    end
end

%%
close all

figure(3)
subplot(221),imagesc(hiv_sc_train_12wk),title('hiv_12wk sc train')
subplot(222),imagesc(hiv_sc_test_12wk),title('hiv_12wk sc test')
subplot(223),imagesc(hiv_fc_train_12wk),title('hiv_12wk fc train')
subplot(224),imagesc(hiv_fc_test_12wk),title('hiv_12wk fc test')

figure(1)
subplot(221),imagesc(hiv_sc_train),title('hiv sc train')
subplot(222),imagesc(hiv_sc_test),title('hiv sc test')
subplot(223),imagesc(hiv_fc_train),title('hiv fc train')
subplot(224),imagesc(hiv_fc_test),title('hiv fc test')


figure(2)
subplot(221),imagesc(hc_sc_train),title('hc sc train')
subplot(222),imagesc(hc_sc_test),title('hc sc test')
subplot(223),imagesc(hc_fc_train),title('hc fc train')
subplot(224),imagesc(hc_fc_test),title('hc fc test')

%%

save('train_test_dataset_agematch_destrieux_1004.mat', 'hiv_sc_train','hiv_sc_test','hiv_fc_train','hiv_fc_test',...
    'hc_sc_train','hc_sc_test','hc_fc_train','hc_fc_test','hiv_sc_train_12wk','hiv_sc_test_12wk','hiv_fc_train_12wk','hiv_fc_test_12wk')

%save('train_test_dataset_agematch_0927_all.mat')
