%% Mean RP amplitude for a specific time window
% ======================================================================= %
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on July 2018. Modified July 2019.
% EXPERIMENT: Timelimit_2018


%{ 

    SCOPE: Script for extracting mean RP amplitude for plotting and stats.

    OUTPUT: datamatrix_premov_20{i,k}, mean_premov_amp_ch20(i,k),
          median_premov_amp_ch20(i, k).
    
    HOW:
    
        1)Load Timelimit_*_subj**_EEG_clean_concat_rej_interp.mat.
        2)Create a new matrix for subjects (11) x conditions (5) and extraxt amp.
        3)Do mean and median of datamatrix_premov_**{i,k}.avg.
        4)Save datamatrix_premov_**, mean_premov_amp_ch**,median_premov_amp_ch**.
        5)Plot mean_all, median_all.

%}

% FIX ME: the loop to load each preprocessed file still not working

%=========================================================================%
%% START of the script 

%% Housekeeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear workspace (if needed)
if input('clear all?  (1/0) ... ')
    clearvars; close all;
end

% set paths (if needed)
BT_setpath

% choose subj & go to the right folder
BT_getsubj

clear LevelAnalysis name numlines prompt subj_folders 

%% More specific paths (maybe set this in the start script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

behavioral_folder= [results_Path, '/Behaviour']; % it can be also current_subj_folder
if ~exist(fullfile(behavioral_folder)); mkdir(fullfile(behavioral_folder)); end;

timeseries_folder= [results_Path, '/Timeseries']; % it can be also current_subj_folder
if ~exist(fullfile(timeseries_folder)); mkdir(fullfile(timeseries_folder)); end;
cd(timeseries_folder);

powerspectra_folder= [results_Path, '/Powerspect']; % it can be also current_subj_folder
    if ~exist(fullfile(powerspectra_folder)); mkdir(fullfile(powerspectra_folder)); end;
    cd(powerspectra_folder);
    
correlation_folder= [results_Path, '/Correlations']; % it can be also current_subj_folder
    if ~exist(fullfile(correlation_folder)); mkdir(fullfile(correlation_folder)); end;  
    
regression_folder= [results_Path, '/Regressions']; % it can be also current_subj_folder
    if ~exist(fullfile(regression_folder)); mkdir(fullfile(regression_folder)); end;
    
statistics_folder= [results_Path, '/Statistics']; % it can be also current_subj_folder
    if ~exist(fullfile(statistics_folder)); mkdir(fullfile(statistics_folder)); end;
      
%% Load timelocked averages from each subject's folders

% Old path:
% 'subj%02d_noHPI/Timeseries/Amplitudes/Baseline_Haggard/timeseries_EEG.mat'
% & 'subj%02d_noHPI/Timeseries/Amplitudes/No_Baseline/timeseries0_EEG.mat';
% 'subj%02d_RP_avg.mat'
 

% Consider only GOOD SUBJECTS 

% With clear Readiness Potential(= negative slope vs flat slope or positive slope).
good_subjs = [2 3 5 6 7 8 10 11 13 15 17 18 20]; % % removed: subj04,subj14,subj19; remove subj 12 & 16 for channel 28 and 30.
    nGoodSubjs = length(good_subjs);
    
OKsubjs= [3 6 7 8 10 13 15 17 18 19 20 21]; % done with the 1µV method assessment 

% like this
% avg2matrix{i,k}= pickupSub(i).avg{k}; %avg_EEG
% stdmatrix{i,k}= pickupStd(i).across_stdev{k}.avg
        
cd(pwd);
for subi=1:nSubjs;
    
    fname_ER= sprintf('subj%02d_TimeS_cond',subi); % or TimeS_condTrl
    pickupER(subi) = load(fname_ER);
   
end

for subi= 1:nSubjs; 
 
    for k= 1:5
        TimeSmatrix{subi,k}= pickupER(subi).avg_cond{k}; % or avg_condTrl or avg_trl or avg_one or st
    end

end

% or load

%% Create a new matrix for subjects (11) x conditions (5) and extraxt amp
% Parameters here:
latHp= [-1 -0.2]; % - 200ms because fo the Point of non-return finding (Schultze-Kraft 2016); -1s not clear 
latExp= [-2 _0.2]; % same but -2s because it's the maximum waiting time for the shortest condition
latEff= [];
% ROI= {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
% chans= 'EEG030'; % or 'EEG020'or 'EEG028' or ROI


datamatrix_premov= struct('ch20',[],'ch28', [],'ch30',[],'ROI', []); 
mean_premov_amp= struct('ch20',[],'ch28', [],'ch30',[],'ROI', []); 
sem_premov_amp= struct('ch20',[],'ch28', [],'ch30',[]);

for i= 1:nSubjs; % OKsubjs
    
    for k= 1:5
        
        %         avg2matrix{i,k}= pickupSub(i).avg{k}; %avg_EEG
        %         stdmatrix{i,k}= pickupStd(i).across_stdev{k}.avg
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         cfg = [];
        %         cfg.latency = [-1 -.2];
        %         cfg.channel = 'EEG020';
        %         datamatrix_premov.ch20{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
        %         mean_premov_amp.ch20(i, k) = mean(datamatrix_premov.ch20{i,k}.avg);
        % %         sem_premov_amp.ch20(i, k) = sem(mean_premov_amp.ch20(i, k),1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         cfg = [];
        %         cfg.latency = [-1 -.2];
        %         cfg.channel = 'EEG028';
        %         datamatrix_premov.ch28{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
        %         mean_premov_amp.ch28(i, k) = mean(datamatrix_premov.ch28{i,k}.avg);
        % %         sem_premov_amp.ch28(i, k) = sem(mean_premov_amp.ch28(i, k),1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cfg = [];
        cfg.latency = latExp;
        cfg.channel = 'EEG030';
        datamatrix_premov.ch30{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
        mean_premov_amp.ch30(i, k) = mean(datamatrix_premov.ch30{i,k}.avg);
        sem_premov_amp.ch30(i, k) = sem(mean_premov_amp.ch30(i, k),1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         cfg = [];
        %         cfg.latency = [-1 -0];%[-1 -.2]
        %         cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
        %         datamatrix_premov.ROI{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
        %         datamatrix_premov.ROI{i,k}= mean(datamatrix_premov.ROI{i,k}.avg,1);
        %         mean_premov_amp.ROI(i, k) = mean(datamatrix_premov.ROI{i,k});
        %         sem_premov_amp.ROI(i, k) = sem(mean_premov_amp.ROI(i, k),1);
    end
end

%% Save

cd(timeseries_folder);

save datamatrix_premov datamatrix_premov; save mean_premov_amp mean_premov_amp 

%% GRANDAVERAGE of mean amplitudes for specific window of time

%Function rule for Recursive sequence (24 Oct)
a=2; r=2;n=5;
s = a*r.^(0:n-1); 

mean_all= struct('ch20',[],'ch28', [],'ch30',[],'ROI', []);
sem_all= struct('ch20',[],'ch28', [],'ch30',[]),'ROI', [];

    
    mean_all.ch20= mean(mean_premov_amp.ch20(:,:));
    sem_all.ch20= sem(mean_premov_amp.ch20(:,:),1);
    
    mean_all.ch28= mean(mean_premov_amp.ch28(:,:));
    sem_all.ch28= sem(mean_premov_amp.ch28(:,:),1);
    
    
    mean_all.ch30= mean(mean_premov_amp.ch30(OKsubjs,:));
    sem_all.ch30= sem(mean_premov_amp.ch30(OKsubjs,:),1);
    
    mean_all.ROI= mean(mean_premov_amp.ROI(OKsubjs,:));
    sem_all.ROI= sem(mean_premov_amp.ROI(OKsubjs,:),1);

cd(timeseries_folder);    
save GAVG_Amp mean all sem_all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END