%% Grandaverage for Timelimit2018 (ft_timelockgrandaverage)
%=========================================================================%
% This routine computes the grandmean of RP amplitudes across participants 
% with the Fiedltrip function 'ft_timelockgrandaverage'.

% AUTHOR: Bianca Trov?.
% DATE: 3 June 2018. Modified: 27 September.

% Outputs:  datamatrix{i,k}, 11x5 cell, all avg for each subj, each cond.
%           Grand_Avg{k}, 1x5 cell, grandavg for each cond, all subjs.


% HOW:
%{ 

1)Put all subjs timelock files in pickupSub(i), 1x11 struct.
2)'Reshape' in a matrix subjects x conditions: datamatrix{i,k}, 11x5 cell.
3)Use it to do grandaverage across conditions in Grand_Avg{k},1x5 cell.
4)save datamatrix and Grand_Avg.
5)plot Grand_Avg with ft_multiplotER.

%}

% FIX ME: add legend to ft_multiplotER
% NOTE: WHICH FILTERS, WHICH BASELINE??

%=========================================================================%
%% Set path (27.09)
clearvars;
% Define some paths 
if strcmp(computer, 'MACI64')% on my laptop
    script_Path= '/Volumes/USB_DISK/TIMELIMIT_backup/SCRIPTS_ANALYSES/MEEG'; % here you find all the scripts for preprocessing/analysing MEEG data
    data_Path = '/Volumes/USB_DISK/TIMELIMIT_backup/MEEG_fif_files'; % = parent_folder: all the raw data are stored here.
    ft_Path = '/Users/bt_neurospin/matlab/FIELDTRIP/fieldtrip-20170405'; % Fieldtrip tools 
    tool_Path= '/Users/bt_neurospin/matlab'; % other useful functions for matlab (written by Aaron)
end

% add some paths
addpath(genpath(script_Path)); % general pre-proc path
addpath(genpath(tool_Path));
addpath(ft_Path); ft_defaults; % start up fieldtrip [NEW]
addpath([ft_Path '/engine']); % start up FT engines [NEW] NOTE: what is this actually doing? 

cd(data_Path); 
parent_folder= data_Path;
subj_folders = dir(fullfile(parent_folder, 'subj*'))
nSubjects= length(subj_folders);
% Consider only good subjects with clear Readiness Potential (= negative
% slope vs flat slope or positive slope).
good_subjects = [2 3 5 6 7 8 10 11 13 15 17 18 20]; % % removed: subj04,subj14,subj19; remove subj 12 & 16 for channel 28 and 30.
nGoodSubjects = length(good_subjects);
timers = [2 3 6 7 9 13 15 18];
nTimeSubjects= length(timers);
non_timers= [1 4 5 8 10 11 12 14 16 17 19 20];
nNon_TSubjects= length(non_timers);
grp1= [4 8 10 12 14 16];
grp2= [1 2 3 5 6 7 9 11 13 15 17];

%% Load timelock files from all the subjects' folders

for subi=1:nSubjects;
    if subi== 1
%                 fname = sprintf('subj%02d_noHPI/Timeseries/Amplitudes/Baseline_Haggard/timeseries_EEG.mat',subi);
        fname = sprintf('subj%02d_noHPI/Timeseries/Amplitudes/No_Baseline/timeseries_EEG.mat',subi);
    else
%                 fname = sprintf('subj%02d/Timeseries/Amplitudes/Baseline_Haggard/timeseries_EEG.mat',subi); %'subj%02d/avg_EEG.mat'
        fname = sprintf('subj%02d/Timeseries/Amplitudes/No_Baseline/timeseries_EEG.mat',subi);
    end
    if isfile(fname)
        pickupSub(subi) = load(fname)
    end
end

%%
% added November 11
% 
% for subi=1:nSubjects;
%     if subi== 1
%         fname = sprintf('subj%02d_noHPI/Timeseries/Amplitudes/Baseline_Haggard/avg_EEG.mat',subi);
%     else
%         fname = sprintf('subj%02d/Timeseries/Amplitudes/Baseline_Haggard/avg_EEG.mat',subi);
%     end
%     if isfile(fname)
%         pickupAvg(subi) = load(fname);
%     end
% end
%%
for subi=1:nSubjects;
    if subi== 1
        fname = sprintf('subj%02d_noHPI/Timeseries/Amplitudes/Baseline_Haggard/std_EEG.mat',subi);
%         fname = sprintf('subj%02d_noHPI/Timeseries/Amplitudes/No_Baseline/std_EEG.mat',subi);
    else
        fname = sprintf('subj%02d/Timeseries/Amplitudes/Baseline_Haggard/std_EEG.mat',subi);
%         fname = sprintf('subj%02d/Timeseries/Amplitudes/No_Baseline/std_EEG.mat',subi);
    end
    if isfile(fname)
        pickupStd(subi) = load(fname);
    end
end

%% Remove more subjects
bad_subjects = [1 4 9 12 14 16 ]; % % removed: subj04,subj09, subj14; remove subj 12 & 16 for channel 28 and 30.
nBadSubjects = length(bad_subjects);

for subi= 1:nBadSubjects
   pickupSub(bad_subjects(subi)).avg = NaN;
end

%% Create a new matrix for subjects (11)or (17) x conditions (5)
avgmatrix=[]; stdmatrix=[];
for subi= 1:nSubjects; %nGoodSubjects
    for k= 1:5
        avgmatrix{subi,k}= pickupSub(subi).avg{k}; %pickupSub(i).avg_EEG{k}
%         avgmatrix{subi,k}= pickupSub(subi).avg{k}; 
        stdmatrix{subi,k}= pickupStd(subi).across_stdev{k}; 
    end
end

%% Do grandaverage across conditions
Grand_Avg=[]; Grand_Std=[];
for k= 1:5
cfg=[];
cfg.preproc.lpfilter='yes';
cfg.preproc.lpfreq= 2;
    Grand_Avg{k}= ft_timelockgrandaverage([],avgmatrix{:,k});
    Grand_Std{k}= ft_timelockgrandaverage([],stdmatrix{:,k});
end

%% Save stuff in the right folder (27.09)

% Create the folder if it doesn't exist already.
% [status, msg] = mkdir('Results'); %IT'S NOT WORKING PROPERLY

%create folder if it doesn't already exist
[status, msg] = mkdir(parent_folder,'Results/Timeseries/Amplitudes/No_Baseline');
cd(fullfile(parent_folder,'Results/Timeseries/Amplitudes/No_Baseline'));

% save datamatrix datamatrix
% save Grand_Avg Grand_Avg 

save avg2matrix avgmatrix; save stdmatrix stdmatrix;
save Grand_2Avg Grand_Avg; save Grand_Std Grand_Std;

%% load
cd(fullfile(parent_folder,'Results/Timeseries/Amplitudes/No_Baseline'));
load('avgmatrix'); load('stdmatrix'); load('Grand_Avg'); load('Grand_Std');

%% PLOT
% multiplot 
cfg=[];
% cfg.preproc.lpfilter='yes';
% cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
% cfg.baseline = 'yes';
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 1.5;
figure
ft_multiplotER(cfg,Grand_Avg{1},Grand_Avg{2},Grand_Avg{3},Grand_Avg{4},Grand_Avg{5});

%% 
cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 1.5;
% cfg.baseline = 'yes';
% cfg.baselinewindow = [-0.005 0.005];
figure
ft_multiplotER(cfg,Grand_Std{1},Grand_Std{2},Grand_Std{3},Grand_Std{4},Grand_Std{5});
% singleplot (to look at single RPs)

%% Plot to check if it shows RP
conditions_all= [1 2 3 4 5];
n_conditions= length(conditions_all);
channels= [20 28 30]; %20= FC1; 28= C3; 30= CZ (EEG cap 60 electrodes)

for j= 1:numel(channels)
    h(j)=figure;
    for condi= 5
        
        selected_data = Grand_Avg{condi}.avg(channels(j),:); % 1001:1500 samples for statistics
        selected_time = Grand_Avg{condi}.time(:);
        h(j,condi)=figure;
        % subplots?
        %         subplot(3,5,15)
        plot(selected_time, selected_data,'LineWidth',3)
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        title(['Grandaverage Readiness potential from condition ' num2str(condi) ', channel ' num2str(channels(j))])
        %         title(['Cond ' num2str(condi) ', chan' num2str(channels(j))])
        xlim([-3.0 1.0])
        %         ylim([-1e-13 3e-13]) %??
        %         hold on
        % filename= ['Grandavg_RP_allconds', '_Chan_' int2str(channels(j)) '.png'];
        
    end
    filename= ['Grandavg_RP_Inf', '_Chan_' int2str(channels(j)) '.png'];
        saveas(h(j),filename); % you will save 15 figures (5x3)
        hold on
        hold off
end

% topoplot (try for Early RP [-1 -0.2] vs Late RP [-0.2 0])
%% regular plot
conditions_all= [1 2 3 4 5];
n_conditions= length(conditions_all);
channels= [20 28 30]; %20= FC1; 28= C3; 30= CZ (EEG cap 60 electrodes)


for j= 1:numel(channels)
    k(j)=figure;
    for condi= 1: n_conditions;
        
        colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]} % Cell array of colors.
        hold on
        p1=plot(Grand_Avg{condi}.time(:),Grand_Avg{condi}.avg(channels(j),:),'LineWidth',2); %this time use samples on the time axes
        hold on
        title(['Grandaverage RP amplitude for short vs long conditions, ' ' Channel ' int2str(channels(j)) ]);
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        legend('2sec','4sec','8ec','16sec','Inf','Location','SouthEast')% COLORS DO NOT MATCH
        hold on
        filename= ['RP_grandavg_Chan_' int2str(channels(j)) '_allconds', '.png'];
        
    end
    saveas(k(j),filename); %HOW TO SAVE ONLY LAST FIGURE FROM A LOOP??
    hold on
    hold off
end

%% END (for now)
%%%%%%%%%%%%%%%%%%%