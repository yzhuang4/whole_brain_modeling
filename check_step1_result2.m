clear all
close all
clc

%figure_dir = '/gpfs/fs1/home/yzhuang4/Downloads/CBIG-0.7.0-Wang2018_MFMem/stable_projects/fMRI_dynamics/Wang2018_MFMem/result_figures';
figure_dir = '/gpfs/fs2/scratch/yzhuang4/CBIG-0.7.0-Wang2018_MFMem/result_figures'
%% 1) load estimated parameters for multiple randomized initializations.

% load healthy control.

num_iterations = 500;
filelist = dir('Estimated_Parameter_hc_500*.mat');
num_files = length(filelist);

hc.Para_E = zeros([138, num_files]);
hc.rrr_z_max = zeros([1, num_files]);
hc.rrr_z = zeros([num_iterations, num_files]);


for idx_file = 1:num_files
    tmp = load(filelist(idx_file).name);
    hc.Para_E(:, idx_file) = tmp.Para_E;
    hc.rrr_z_max(idx_file) = tmp.rrr_z_max;
    hc.rrr_z(:, idx_file) = tmp.rrr_z;
    hc.num_files = num_files;
end

% load hiv baseline
filelist = dir('Estimated_Parameter_hiv_500*.mat');
num_files = length(filelist);

hiv.Para_E = zeros([138, num_files]);
hiv.rrr_z_max = zeros([1, num_files]);
hiv.rrr_z = zeros([num_iterations, num_files]);

for idx_file = 1:num_files
    tmp = load(filelist(idx_file).name);
    hiv.Para_E(:, idx_file) = tmp.Para_E;
    hiv.rrr_z_max(idx_file) = tmp.rrr_z_max;
    hiv.rrr_z(:, idx_file) = tmp.rrr_z;
    hiv.num_files = num_files;

end


filelist = dir('Estimated_Parameter_hiv_12wk_500*.mat');
num_files = length(filelist);

wk12.Para_E = zeros([138, num_files]);
wk12.rrr_z_max = zeros([1, num_files]);
wk12.rrr_z = zeros([num_iterations, num_files]);

for idx_file = 1:num_files
    tmp = load(filelist(idx_file).name);
    wk12.Para_E(:, idx_file) = tmp.Para_E;
    wk12.rrr_z_max(idx_file) = tmp.rrr_z_max;
    wk12.rrr_z(:, idx_file) = tmp.rrr_z;
    wk12.num_files = num_files;

end

%%  2) step1 result,shows simulation results perform better in hiv bsl > wk12 > hc.

clf
figure(1),

hold on

for i = 1:10 %size(hc.rrr_z,2)
    f1 = scatter(1:num_iterations, hc.rrr_z(:,i), ...
        5, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerFaceColor','r');
    hold on
end
hold on

for i = 1:10 %size(wk12.rrr_z,2)
    f3 = scatter(1:num_iterations, wk12.rrr_z(:,i), ...
        5, 'filled','MarkerFaceAlpha', 0.5, 'MarkerFaceColor','g');
    hold on
end

for i = 1:10 %size(hiv.rrr_z,2)
    f2 = scatter(1:num_iterations, hiv.rrr_z(:,i), ...
        5, 'filled','MarkerFaceAlpha', 0.5,'MarkerFaceColor', 'b');
    hold on
end

hc_mean_rrrz = mean(hc.rrr_z,2);
hiv_mean_rrrz = mean(hiv.rrr_z,2);
wk12_mean_rrrz = mean(wk12.rrr_z,2);
%plot(hc_mean_rrrz,'Color','r','LineWidth',2)
%plot(hiv_mean_rrrz,'Color','b','LineWidth',2)
%plot(wk12_mean_rrrz,'Color','g','LineWidth',2)


xlabel('number of iterations')
xlim([-10,num_iterations + 50])
ylim([0.001,max(max(hiv.rrr_z))+0.1])
ylabel('similarity')
title('similarity z-score between simulated and empirical FC')
legend([f1 f2 f3], 'HC','HIV+BSL', 'HIV+12wsk', 'Location', 'southeast')


hold off

%saveas(figure(1), [figure_dir, '/agematch_relax-mfm_result_500iters_rmltnp.fig'])
%saveas(figure(1), [figure_dir, '/agematch_relax-mfm_result_500iters_rmltnp.pdf'])


%%

%% 

figure(2)
f1 = scatter(1:hc.num_files, hc.rrr_z_max,150, 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor','r');
hold on
plot(1:hc.num_files, hc.rrr_z_max,'r')

plot(1:hiv.num_files, hiv.rrr_z_max,'b')
f2 = scatter(1:hiv.num_files, hiv.rrr_z_max,150, 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor','b');
plot(1:wk12.num_files, wk12.rrr_z_max,'g')
f3 = scatter(1:wk12.num_files, wk12.rrr_z_max,150, 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor','g');
ylabel('maximum similartiy')
xlabel('index of random initializations')
legend([f1 f2 f3], 'HC','HIV+BSL','HIV+12wk', 'Location', 'southeast')
hold off
%saveas(figure(2), [figure_dir, '/agematch_relax-mfm_result_500iters2_rmltnp.fig'])
%saveas(figure(2), [figure_dir, '/agematch_relax-mfm_result_500iters2_rmltnp.pdf'])


figure(3)
rrr_z_max_3grps = [hiv.rrr_z_max, wk12.rrr_z_max, hc.rrr_z_max];
grp = [1*ones(1,hiv.num_files), 2*ones(1,wk12.num_files), 3*ones(1,hc.num_files)];
%grp = [repmat({'hiv bsl'},1,hiv.num_files),repmat({'hiv 12wk'},1,wk12.num_files),repmat({'hc'},1,hc.num_files)];

h = boxplot(rrr_z_max_3grps, grp, 'Labels',{'HIV BSL','HIV Wk12','HC'}, 'Colors','k')
hold on
set(h, 'linewidth', 2)
scatter(1*ones(1,hiv.num_files),hiv.rrr_z_max, 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor','b')
scatter(2*ones(1,wk12.num_files), wk12.rrr_z_max,'filled', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor','r')
scatter(3*ones(1,hc.num_files), hc.rrr_z_max,'filled', 'MarkerFaceAlpha', 0.4, 'MarkerFaceColor','g')

ylabel('maximum similartiy')
xlabel('group')
hold off
%saveas(figure(3), [figure_dir, '/agematch_relax-mfm_result_500iters3_rmltnp.fig'])
%saveas(figure(3), [figure_dir, '/agematch_relax-mfm_result_500iters3_rmltnp.pdf'])


%% 3) use the top 5 iters which have the best rrr_z_max, average the Para_E, for next step.

[hc.rrr_z_max_top5, hc.rrr_z_max_top5_idx] = maxk(hc.rrr_z_max, 5);
[hiv.rrr_z_max_top5, hiv.rrr_z_max_top5_idx] = maxk(hiv.rrr_z_max, 5);
[wk12.rrr_z_max_top5, wk12.rrr_z_max_top5_idx] = maxk(wk12.rrr_z_max, 5);

hc.Para_E_avg = mean(hc.Para_E(:,hc.rrr_z_max_top5_idx),2);
hiv.Para_E_avg = mean(hiv.Para_E(:,hiv.rrr_z_max_top5_idx),2);
wk12.Para_E_avg = mean(wk12.Para_E(:,wk12.rrr_z_max_top5_idx),2);






% there is no diffrence in regard of Para_E_avg, rrr_z_max between 0717 and
% 0708. So step2 results using 0708 is valid.

%save('500iters_final_estimate_Paras_agematched_1003.mat','hc','hiv','wk12');
