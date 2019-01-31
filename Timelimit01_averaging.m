%% AVERAGING ERPs/ERFs for Timelimit2018 (ft_timelockanalysis)
% ======================================================================= %
% Script for post-processing EEG/MEG data from the Timelimit project(2018).   
% This routine performs timelockanalysis on pre-processed data. It takes
% for example Timelimit_*_subj**_EEG_clean_concat_rej_interp.mat then...
%
% AUTHOR: Bianca Trovo (PhD student)
% CREATED: 25th May 2018. Modified: 28 May, 29 May,30 May.

% Output: avg_EEG{i}, 1x5 cell.

% HOW:

%{ 

1)Load data already preprocessed with myft_preProcess_and_clean_data.m
    (Aaron's script) and saved as
    Timelimit_*_subj**_EEG_clean_concat_rej_interp.mat.
2)Apply baseline?
3)Do ft_timelockanalysis.
4)Save avg_EEG.
5)Plot.

%}

% FIX ME: the loop to load each preprocessed file still not working
% NOTE: 

%=========================================================================%
%% BEGIN: set paths (27.09)
fprintf('\n BEGIN: clear variables and set correct path \n')
clearvars;

% Define some paths 
if strcmp(computer, 'MACI64')% on my laptop
    script_Path= '/Volumes/USB DISK/TIMELIMIT_backup/SCRIPTS_ANALYSES/MEEG'; % here you find all the scripts for preprocessing/analysing MEEG data
    data_Path = '/Volumes/USB DISK/TIMELIMIT_backup/MEEG_fif_files'; % = parent_folder: all the raw data are stored here.
    ft_Path = '//Users/bt_neurospin/matlab/FIELDTRIP/fieldtrip-20170405'; % Fieldtrip tools 
    tool_Path= '/Users/bt_neurospin/matlab'; % other useful functions for matlab (written by Aaron)
end

% add some paths
addpath(genpath(script_Path)); % general pre-proc path
addpath(genpath(tool_Path));
addpath(ft_Path); ft_defaults; % start up fieldtrip [NEW]
addpath([ft_Path '/engine']); % start up FT engines [NEW] NOTE: what is this actually doing? 

% Move to the right folder/path

parent_folder= data_Path; 
% subj_folders = dir(fullfile(parent_folder, 'subj*'));
% current_subj_folder = fullfile(parent_folder, subj_folders(subjnum).name);
cd(parent_folder);

%% Load data from current subject/folder

subj_list= [3 4 5 6 7 9 10 11 12 13 14 15 16];
% nSubjs= length(subj_list);
% 
% for i = 1:nSubjs;
%     filename = sprintf('Time_Limit_2_subj%02d_EEG_clean_concat_rej_interp.mat',subj_list(i));
%     load(fullfile(datapath, filename))
% end
% load([current_subj_folder,'TimeLimit_2_subj*_EEG_clean_concat_rej_interp/.mat'])

subj_folders = dir(fullfile(parent_folder, 'subj*'))
nSubjects= length(subj_folders);
% resultsfolder = fullfile(parent_folder,'/Results');
% cd(resultsfolder);

% for subj=2:nSubjects;
%     current_subject = fullfile(parent_folder, subj_folders(subj).name)
%     load([current_subject,'/TimeLimit_2_subj*_EEG_clean_concat_rej_interp.mat'])
% end DOESNT WORK

for i = 2:nSubjects;
    filename = sprintf('Time_Limit_2_subj%02d_EEG_clean_concat_rej_interp.mat',subj_folders.name(i));
    load(fullfile(parent_folder, filename))
end %DOESN'T WORK


%% Load files
% you should loop it across subjects in each folder. Try first with one
% subject at the time.
% FIX-ME
load TimeLimit_2_subj08_EEG_clean_concat_rej_interp.mat
%% Re-organize data, preprocess and apply baseline

% datapath: '/Volumes/USB DISK/TIMELIMIT_backup/MEEG_fif_files'
% subj_list= [3 4 5 6 7 9 10 11 12 13 14 15 16];
% nSubjs= length(subj_list);
% 
% for i = 1:nSubjs;
%     filename = sprintf('Time_Limit_2_subj%02d_EEG_clean_concat_rej_interp.mat',subj_list(i));
%     load(fullfile(datapath, filename))
% end

% cond= [TRIALS.cond]; %we put all the conditions in a row
% un_conds = unique(cond); %we take the number of conditions with no repetition
 
%% LOAD data cleaned from cognitive artifacts

load('BADcleanedtrls.mat');
good_trls = setdiff([1:length(DATA_REJ_INTERP.trial)],allbad_trls); % length of the number of trials kept 


cond= [TRIALS.cond]; %we put all the conditions in a row
cond(cond==32) = Inf;
un_conds = unique(cond);
newcond= cond(good_trls);

cfg= [];
cfg.demean = 'yes';
dataprep= ft_preprocessing(cfg,DATA_REJ_INTERP);
% cfg.baselinewindow = [-0.05 0.05]; % to get rid of movement stuff.Following Aaron/Patrick's suggestion. Better not to use it now. 
% dataprep_bl = ft_preprocessing(cfg,DATA_REJ_INTERP);

%% Do ERPs/ERFs average for each condition


for i = 1:length(un_conds)
cfg=[];
cfg.trials= find(cond == un_conds(i));
avg_EEG{i}= ft_timelockanalysis(cfg,dataprep); % dataprep_bl if you apply baseline (see above)
end

save avg_EEG avg_EEG

%% Plot the data
% COLOR CODE
% Blue: 2 sec [trigger: 2]
% Red: 4 sec [trgg: 4]
% Green: 8 sec [trgg: 8]
% Black: 16 sec [ trgg: 16]
% Yellow: Inf [trgg: 32]

cfg= [];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.preproc.lpfilter= 'yes';
cfg.preproc.lpfreq= 20;
% cfg.title= 
% cfg.linewidth= 2.0;
ft_multiplotER(cfg,avg_EEG{1},avg_EEG{2},avg_EEG{3},avg_EEG{4},avg_EEG{5});

cfg= [];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.preproc.lpfilter= 'yes';
cfg.preproc.lpfreq= 35;
cfg.channel= [20]
cfg.linewidth= 2.0;
ft_singleplotER(cfg,avg_EEG{1},avg_EEG{2},avg_EEG{3},avg_EEG{4},avg_EEG{5});

% only one condition plotted 
cfg= [];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.preproc.lpfilter= 'yes';
cfg.preproc.lpfreq= 35;
% cfg.linewidth= 2.0;
ft_multiplotER(cfg,avg_EEG{5}); %for the Readiness Potential

cfg= [];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.preproc.lpfilter= 'yes';
cfg.preproc.lpfreq= 35;
cfg.channel= [27,28,29,30]
cfg.linewidth= 2.0;
ft_singleplotER(cfg,avg_EEG{5});

%% END of the script

fprintf('\n END: DATA SAVED \n')