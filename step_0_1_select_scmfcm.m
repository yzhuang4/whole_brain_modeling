%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step0_1_select_scmfcm.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Age matched group, remove some severe headmotion subjects, rerun after
%slice timing correction

clear all
close all
clc
%% load sc and fc
cd('/Users/yuchuan/Box Sync/Connectome/')
load scm.mat
load fcm_st.mat
load demographic.mat
age = [demo_subjID{:,2}];
%% get the exclusion list, 1 - keep, 0 - exclude
filename = '/Users/yuchuan/Box Sync/Connectome/scripts/agematched_rerun/age_match_exclude.txt';
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
        scm(:,:,count,:,:) = connectome_desikan_2(:,:,i,:,:);
        fcm(:,:,count,:) = fcm_desikan_2(:,:,i,:);
        zfcm(:,:,count,:) = fcm_z_desikan_2(:,:,i,:);
        count = count + 1;
    end
end
%%
num_hiv = 46;
num_hc = 46;
num_subject = 92;
save('age_matched_data_st.mat','scm','fcm','zfcm','subj_list_agematch','num_hiv','num_hc','num_subject')

%% load subject list from 'train_test_dataset_agematch_0524_all.mat'
load('train_test_dataset_agematch_0524_all.mat', 'train_hc','train_hiv','train_wk12','test_hc','test_hiv','test_wk12');
%% Use the index to calculate the mean SC_train, mean FC_train, mean SC_test, FC_test
% only extract the SC including the cortical regions 

num_rois = 68;

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
subplot(223),imagesc(hiv_fc_train_12wk),caxis([0 1]),title('hiv_12wk fc train')
subplot(224),imagesc(hiv_fc_test_12wk),caxis([0 1]),title('hiv_12wk fc test')

figure(1)
subplot(221),imagesc(hiv_sc_train),title('hiv sc train')
subplot(222),imagesc(hiv_sc_test),title('hiv sc test')
subplot(223),imagesc(hiv_fc_train),caxis([0 1]),title('hiv fc train')
subplot(224),imagesc(hiv_fc_test),caxis([0 1]),title('hiv fc test')


figure(2)
subplot(221),imagesc(hc_sc_train),title('hc sc train')
subplot(222),imagesc(hc_sc_test),title('hc sc test')
subplot(223),imagesc(hc_fc_train),caxis([0 1]),title('hc fc train')
subplot(224),imagesc(hc_fc_test),caxis([0 1]),title('hc fc test')
%%

save('train_test_dataset_agematch_0625.mat', 'hiv_sc_train','hiv_sc_test','hiv_fc_train','hiv_fc_test',...
    'hc_sc_train','hc_sc_test','hc_fc_train','hc_fc_test','hiv_sc_train_12wk','hiv_sc_test_12wk','hiv_fc_train_12wk','hiv_fc_test_12wk')
