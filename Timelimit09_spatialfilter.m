%% Spatial Filter script
%=========================================================================%
% AUTHORs: Aaron & Yvonne; Bianca modified
% DATE: February/March 2019

%{

    OUTPUT that you need: SF_timecourses_bl (coming from trl + baseline),
    msf and d.


%}
%=========================================================================%
%% START of the script 
%% Set path (27.09) (change it to a main script?)
clearvars;
% Define some paths 
if strcmp(computer, 'MACI64')% on my laptop
%     script_Path= '/Volumes/USB_DISK/TIMELIMIT_backup/SCRIPTS_ANALYSES/MEEG'; % here you find all the scripts for preprocessing/analysing MEEG data
    data_Path = '/Volumes/USB_DISK/TIMELIMIT_backup/MEEG_fif_files'; % = parent_folder: all the raw data are stored here.
    ft_Path = '/Users/bt_neurospin/matlab/FIELDTRIP/fieldtrip-20170405'; % Fieldtrip tools 
    tool_Path= '/Users/bt_neurospin/Repos'; % where all the scripts should be
end

% add some paths
% addpath(genpath(script_Path)); % general pre-proc path
addpath(genpath(tool_Path));
addpath(ft_Path); ft_defaults; % start up fieldtrip [NEW]
addpath([ft_Path '/engine']); % start up FT engines [NEW] NOTE: what is this actually doing? 

%% Choose the subject
prompt={'Which participant do you want to look at?'};
name='Subject number';
numlines=1;
answer=inputdlg(prompt,name,numlines);
subjnum= str2double(answer);

%% Move to the right folder/path

parent_folder=  data_Path;
subj_folders = dir(fullfile(parent_folder, 'subj*')); % 'Subject'
current_subj_folder = fullfile(parent_folder, subj_folders(subjnum).name);
cd(current_subj_folder);

%% Load preprocessed files

if subjnum== 1 || subjnum== 18 || subjnum== 19 || subjnum== 20 || subjnum== 21 || subjnum== 22
    load(sprintf('TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp',subjnum))
else
    load(sprintf('TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp',subjnum))
end


%%
% if subjnum== 17
%     cfg.channel= [1:17 19:32 34:60]
%     DATA_REJ_INTERP= ft_preprocessing(cfg,DATA_REJ_INTERP)
% end

%% Clean data more (need a function?)

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

% just to visualize and find the peak 
gavg = ft_timelockanalysis([],DATA_REJ_INTERP);
cfg.layout = 'eeg_64_NM20884N.lay';
figure
ft_multiplotER(cfg,gavg);

% look at the data before applying SF
TrialAvg= mean(DATA_CUBE,3);
elec_clust= [20,21,29,30,31,39,40]; % cluster of 7 channels selected around Cz
clustavg= mean(TrialAvg(elec_clust,:));

%% Plot to see if there is a nice RP
tAx = [0:2000]./SR - 3;
figure; h=plot(tAx, clustavg)
xlabel('Time (s)')
ylabel('mean Amplitude (\muV)')
title(['Subj ' num2str(subjnum) ', avg channels 20,21,29,30,31,39,40'])
legend('All conditions', 'Location','NorthWest')

%  save figure for further comparisons
filename= [sprintf('Readiness_Potential_subj%02d', subjnum) '.png'];
cd(parent_folder)
saveas(h,filename)

%% ERROR

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

timeseries_folder= [current_subj_folder,'/Timeseries'];
if ~exist([timeseries_folder, '/SpatialFilter']); mkdir([timeseries_folder, '/SpatialFilter']); end
cd(fullfile(timeseries_folder,'/SpatialFilter'));

save SF_results_newBl Ind* tWin nTrials trl msf d blRange SF_timecourses_bl tAx subjnum

%% Tell me what I have done so far
disp(['subj ' num2str(subjnum) ' done'])

close all

%% End of the script