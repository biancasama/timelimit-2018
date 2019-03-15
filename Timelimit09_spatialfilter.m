%% Compute Spatial Filter within subject
%=========================================================================%
% AUTHORs: Aaron & Yvonne; Bianca modified
% DATE: February/March 2019

%{

    OUTPUT that you need: SF_timecourses_bl (coming from trl + baseline),
    msf and d.


%}
%=========================================================================%
%% START of the script 

%% Housekeeping

% clear workspace (if needed)
if input('clear all?  (1/0) ... ')
    clearvars; close all;
end

% set paths (if needed)
BT_setpath

% choose subj & go to the right folder
BT_getsubj

%% Load preprocessed files

if subjnum== 1 || subjnum== 18 || subjnum== 19 || subjnum== 20 || subjnum== 21 || subjnum== 22
    load(sprintf('TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp',subjnum))
else
    load(sprintf('TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp',subjnum))
end

%% Clean data more (need a function?)

cfg=[];
cfg.trials = good_trls;
% cfg.demean = 'yes';
cfg.preproc.lpfilter='yes';
cfg.preproc.lpfreq= 50; %notch filter
DATA_REJ_INTERP= ft_preprocessing(cfg,DATA_REJ_INTERP);

% per condition 
DATA_cond=[];
for condi = 1:length(un_conds)
    cfg=[];
    cfg.trials= find(newcond == un_conds(condi));
    DATA_cond{condi} = ft_preprocessing(cfg,DATA_REJ_INTERP);
%     avg{condi} = ft_timelockanalysis(cfg,DATA_bl);
%     mean_avg{condi} = avg{condi}.avg;
end

%% Do Spatial Filter here

% define parameters here (sampling Rate, Timelock, indexes condition)
SR=500;
t0= 3;
Ind2 = [TRIALS.cond]' == 2;
Ind4 = [TRIALS.cond]' == 4;
Ind8 = [TRIALS.cond]' == 8;
Ind16 = [TRIALS.cond]' == 16;
Ind32 = [TRIALS.cond]' == 32;


% HERE you make the data in 3D structure
DATA_CUBE = myft_ftStruct2dataCube(DATA_REJ_INTERP);


DATA_CUBEcond=[];
for condi = 1:length(un_conds)
    
    DATA_CUBEcond{condi} = myft_ftStruct2dataCube(DATA_cond{condi});
    
end

% just to visualize and find the peak 
gavg = ft_timelockanalysis([],DATA_REJ_INTERP);
cfg.layout = 'eeg_64_NM20884N.lay';
% % cfg=[];
% cfg.preproc.lpfilter='yes';
% cfg.preproc.lpfreq= 50;
figure
ft_multiplotER(cfg,gavg);

% look at the data before applying SF
TrialAvg= mean(DATA_CUBE,3);
elec_clust= [20,21,29,30,31,39,40]; % cluster of 7 channels selected around Cz
clustavg= mean(TrialAvg(elec_clust,:));

TrialAvg= mean(DATA_CUBEcond{5},3);
elec_clust= [20,21,29,30,31,39,40]; % cluster of 7 channels selected around Cz
clustavg= mean(TrialAvg(elec_clust,:));

%% Plot to see if there is a nice RP
tAx = [0:2000]./SR - 3;
figure; H=plot(tAx, clustavg);
xlabel('Time (s)');
ylabel('mean Amplitude (\muV)');
title(['Subj ' num2str(subjnum) ', avg channels 20,21,29,30,31,39,40']);
% legend('All conditions', 'Location','NorthWest');
legend('Infinity condition', 'Location','NorthWest');

%  save figure for further comparisons
filename= [sprintf('Readiness_Potential_subj%02d', subjnum) '.png'];
% cd(parent_folder)
cd(results_Path);
saveas(H,filename);

%  save figure for further comparisons
filename= [sprintf('Readiness_Potential_Inf_subj%02d', subjnum) '.png'];
% cd(parent_folder)
cd(results_Path);
saveas(H,filename);

%% ERROR TOPOPLOT

% cfg = [];                            
% % cfg.xlim = [0.3 0.5];                
% % cfg.zlim = [0 6e-14];                
% cfg.layout = 'eeg_64_NM20884N.lay';            
% cfg.parameter = 'individual'; % the default 'avg' is not present in the data
% figure; ft_topoplotER(cfg,clustavg'); colorbar

% cfg.layout = 'eeg_64_NM20884N.lay';
% figure; topoplot(cfg,clustavg');
% title(['Subj ' num2str(subjnum) ', Topography RP avg channels 20,21,29,30,31,39,40']);
% filename3= [sprintf('Readiness_Potential_topography_subj%02d', subjnum) '.png'];
% cd(parent_folder)
% saveas(k,filename3)

%% Define window for Spatial Filter based on the peak amplitude

Peak_Lat= cursor_info.Position(1);
Peak_Win= [Peak_Lat-0.05 Peak_Lat+0]; %[Peak_Lat-0.05 Peak_Lat+0.05];

% use the time difference to get the spatial filter
tWin= {Peak_Win-1, Peak_Win}; % define tWin based on grand average for that subject

%% Quantitavely assess if there is a RP (microVolt decrease method)
% Alternative method to Slope analysis

M1=mean(clustavg(seconds2samples(tWin{1},t0,SR)));
M2=mean(clustavg(seconds2samples(tWin{2},t0,SR)));
if M2-M1 < -1e-6 % 1 microVolt
    display('OK')
else
    display('NOT OK')
end

%% SPATIAL FILTER here

nTrials = size(DATA_CUBE,3);
IndAll = true(nTrials,1);
[trl,msf] = ems_ncond(DATA_CUBE,IndAll,[], @spf_basic_2time_diff,tWin,t0,SR);

%% using the peak to get the spatial filter
% d = mean(mean(DATA_CUBE(:,seconds2samples(tWin{2},t0,SR),:),3),2);
% d = d ./ norm(d);
% for i=1:162
%     RESULT(i,:) = -(DATA_CUBE(:,:,i)' * d)';
% end
% plot(tAx,mean(RESULT))
% cfg.layout = 'eeg_64_NM20884N.lay';
% figure;topoplot(cfg,d);
% colorbar

% blRange = [0 0.050]; % RP peak
blRange= Peak_Win;
% blRange = [0.225 0.275]; % somato-motor potential
SF_timecourses_bl = baseline_correct(trl,SR,3,blRange);

%% Plot the time course after the Spatial Filter
figure; k=plot(tAx,mean(SF_timecourses_bl))
xlabel('Time (s)')
ylabel('mean Amplitude (\muV)')
title(['Subj ' num2str(subjnum) ', Spatial filter, baseline [ ' num2str(Peak_Win) ' ]'])
legend('All conditions', 'Location','NorthWest')

% save figure for further comparisons
filename= [sprintf('Spatial_filter_subj%02d', subjnum) '.png'];
cd(parent_folder)
saveas(k,filename)

%% Plot the topography after the Spatial Filter

meanSF = mean(msf,2);
cfg.layout = 'eeg_64_NM20884N.lay';
j=figure;topoplot(cfg,meanSF);
title(['Subj ' num2str(subjnum) ', Topography Spatial filter, baseline [ ' num2str(Peak_Win) ' ]']);
colorbar

% save figure for further comparisons
filename= [sprintf('Spatial_filter_topography_subj%02d', subjnum) '.png'];
cd(parent_folder)
saveas(j,filename)

%% Save results for each participants in proper folder 

if input('Save SPATIAL FILTER results? ')
    timeseries_folder= [current_subj_folder,'/Timeseries'];
    if ~exist([timeseries_folder, '/SpatialFilter']); mkdir([timeseries_folder, '/SpatialFilter']); end
    cd(fullfile(timeseries_folder,'/SpatialFilter'));

    save SF_results_newBl Ind* tWin nTrials trl msf d blRange SF_timecourses_bl tAx subjnum
end
%% Tell me what I have done so far
disp(['subj ' num2str(subjnum) ' done'])

close all

%% End of the script