%% CLUSTER-based permutation tests on ERPs
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on April/May 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE:  

    OUTPUT: 

    FIXME: 

%}
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
    
regression_folder2= [regression_folder, '/RegressionsContr']; % it can be also current_subj_folder
    if ~exist(fullfile(regression_folder2)); mkdir(fullfile(regression_folder2)); end;
    
statistics_folder= [results_Path, '/Statistics']; % it can be also current_subj_folder
    if ~exist(fullfile(statistics_folder)); mkdir(fullfile(statistics_folder)); end;
    
alpha_folder= [statistics_folder, '/alpha']; % it can be also current_subj_folder
    if ~exist(fullfile(alpha_folder)); mkdir(fullfile(alpha_folder)); end;
    
prob_folder= [statistics_folder, '/prob']; % it can be also current_subj_folder
    if ~exist(fullfile(prob_folder)); mkdir(fullfile(prob_folder)); end;

%% Load

% Correlations

% load 'R_all_ER1'; load 'R_all_ER2'; 
% load 'R_all_ER1Log'; load 'R_all_ER2Log';
load 'R_all_ER1Log_semiLogX'; load 'R_all_ER2Log_semiLogX'; % correggi
load 'R_all_ER1Contr_semiLogX'; load 'R_all_ER2Contr_semiLogX';

% Regressions
cd(regression_folder); cd(regression_folder2);
% load 'B_all_ER1'; load 'B_all_ER2'; 
% load 'B_all_ER1Log'; load 'B_all_ER2Log';
load 'B_all_ER1Log_semiLogX'; load 'B_all_ER2Log_semiLogX';

load 'B_all_ER1Contr_semiLogX'; load 'B_all_ER2Contr_semiLogX';

%% Parameters

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s

Subjs= 1:22;
RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
BadRPs= Subjs(~ismember(1:numel(Subjs),RPsubjs));
NonTimers= [1 4 5 8 10 16 22];
Timers= Subjs(~ismember(1:numel(Subjs),NonTimers));

%% Randomly choose some participants to exclude from analysis to equalize groups in design matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GoodRP vs Timers
clear Rand;

Rand = Timers(end,randperm(size(Timers,2), 3));
excl= sort(Rand);
Timesubjs= Timers(~ismember(Timers,excl));
% check
size(Timesubjs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timers vs NonTimers
clear Rand;

Rand = Timers(end,randperm(size(Timers,2), 8));
excl= sort(Rand);
Timesubjs= Timers(~ismember(Timers,excl));

% check
size(Timesubjs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GoodRP vs BadRP
clear Rand;

Rand = RPsubjs(end,randperm(size(RPsubjs,2), 2));
excl= sort(Rand);
RPsubjs= RPsubjs(~ismember(RPsubjs,excl));

% check
size(RPsubjs);

%% 1) Time-series CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nSubjs= length(RPsubjs);

% R_all=R_all';

R_some= {};
for subi=1:nSubjs;
    R_some{subi}= R_all{RPsubjs(subi)}
end

% R_some= R_some';

R_some= {};
for subi=1:nSubjs;
    R_some{subi}= R_all{RPsubjs(subi)}
end

isequal(R_all{13},R_some{6}); % IT WORKS but you should automatize it

%% Run cluster test for CORRELATIONS

                % ROI for RP
%                 chans = {'EEG20','EEG21','EEG29','EEG30','EEG31','EEG39','EEG40'};
                % single-channel
                % chans= {'EEG030'}; 
                chans= {'EEG*'};
                latency= [-2  -0.2]; %-200ms = point of no return;
                nSubjs = 12;
                % Prepare the 'design'. Since Fieldrip does not know the the simple t-test
                % against 0, we have to use the paired t-test (against 0).
                design = zeros(2, 2*nSubjs);
                design(1, 1:nSubjs)= 1;               % condition 1 (data)
                design(1, (nSubjs+1):(2*nSubjs))= 2;        % condition 2 (0, value against which we compare the data).
                design(2, 1:nSubjs) = 1:nSubjs;          % Indicate the subject number (for pairing, alghough not useful here)
                design(2, (nSubjs+1):(2*nSubjs)) = 1:nSubjs;

                % Make *fake data*
                R_contr= {};
                for subi=1:nSubjs;

                    templateER.avg = R_some{subi}.avg; % some or all
                    R_contr{subi}= templateER;
                    R_contr{subi}.avg= zeros(nChans,nTimes);
                end

                cfg = [];
                cfg.method = 'montecarlo';
                cfg.channel = chans;
                cfg.latency= latency;
                cfg.avgoverchan = 'yes';
                % cfg.avgovertime= 'yes';
                cfg.correctm = 'cluster';
                cfg.clusteralpha = 0.05; % lower to 0.01
                cfg.correcttail = 'alpha'; %only for two sided t-tests
                cfg.statistic = 'ft_statfun_depsamplesT';
                % cfg.statistic = 'ft_statfun_indepsamplesF'; % if you want
                % to run the ANOVA
                cfg.design = design;
                cfg.ivar = 1;
                cfg.uvar = 2;
                cfg.numrandomization = 500; % to try put 500 and final put 1000

                % elie's way
                cfg_neighb=[];
                % cfg_neighb.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
                cfg_neighb.layout= 'eeg_64_NM20884N.lay';
                cfg_neighb.method = 'distance';
                cfg_neighb.neighbourdist = .18;
                cfg_neighb.sens= R_all{2}.elec;
                neighbs = ft_prepare_neighbours(cfg_neighb);
                cfg.neighbours = neighbs;

                % STATS  
                stat = ft_timelockstatistics(cfg, R_some{:}, R_contr{:}); %R_some or R_all

                cd(statistics_folder);
                save('stat_corr_ER2Contr_semiLog_nr_12', 'stat','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_corr_ER1semiLog_hp_22','stat','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s
                
               % Check if there are significant clusters after correction
                stat.posclusters(1); stat.negclusters;
                
%% 2) Time-series: REGRESSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subgroups

clear nSubjs;
nSubjs= length(RPsubjs);

B_all=B_all';

B_good= {};
for subi=1:nSubjs;
    B_good{subi}= B_all{RPsubjs(subi)}
end

% B_some= B_some';

isequal(B_all{13},B_some{6}); % IT WORKS but you should automatize it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear nSubjs;

nSubjs= length(BadRPs);

B_bad= {};
for subi=1:nSubjs;
    B_bad{subi}= B_all{BadRPs(subi)}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear nSubjs;

nSubjs= length(Timesubjs);

B_time= {};
for subi=1:nSubjs;
    B_time{subi}= B_all{Timesubjs(subi)}
end

isequal(B_good{4},B_time{4});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear nSubjs;

nSubjs= length(NonTimers);

B_NonTime= {};
for subi=1:nSubjs;
    B_NonTime{subi}= B_all{NonTimers(subi)}
end

%% Run cluster test for REGRESSIONS


                % ROI for RP
%                 chans = {'EEG20','EEG21','EEG29','EEG30','EEG31','EEG39','EEG40'};
                % single-channel
                % chans= {'EEG030'}; 
                chans= {'EEG*'};
                latency= [-2  -0.2]; %-200ms = point of no return;
                nSubjs = 22;
                % Prepare the 'design'. Since Fieldrip does not know the the simple t-test
                % against 0, we have to use the paired t-test (against 0).
                design = zeros(2, 2*nSubjs);
                design(1, 1:nSubjs)= 1;               % condition 1 (data)
                design(1, (nSubjs+1):(2*nSubjs))= 2;        % condition 2 (0, value against which we compare the data).
                design(2, 1:nSubjs) = 1:nSubjs;          % Indicate the subject number (for pairing, alghough not useful here)
                design(2, (nSubjs+1):(2*nSubjs)) = 1:nSubjs;

                % Make *fake data*
                B_contr= {};
                for subi=1:nSubjs;

                    templateER.avg = B_all{subi}.avg;
                    B_contr{subi}= templateER;
                    B_contr{subi}.avg= zeros(nChans,nTimes)
                end

                cfg = [];
                cfg.method = 'montecarlo';
                cfg.channel = chans;
                cfg.latency= latency;
%                 cfg.avgoverchan = 'yes';
                % cfg.avgovertime= 'yes';
                cfg.correctm = 'cluster';
                cfg.clusteralpha = 0.05; % 0.05 or 0.01
%                 cfg.alpha = 0.01;
                cfg.correcttail = 'prob'; % IMPORTANT: if you set alpha you should multiply your stat.prob x 2 
                cfg.statistic = 'ft_statfun_depsamplesT';
                % cfg.statistic = 'ft_statfun_indepsamplesF'; % if you want
                % to run the ANOVA
                cfg.design = design;
                cfg.ivar = 1;
                cfg.uvar = 2;
                cfg.numrandomization = 500; % for the paper = 1000.
                
                % elie's way
                cfg_neighb=[];
                % cfg_neighb.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
                cfg_neighb.layout= 'eeg_64_NM20884N.lay';
                cfg_neighb.method = 'distance';
                cfg_neighb.neighbourdist = .18;%.18 %.13
                cfg_neighb.sens= B_all{2}.elec;
                neighbs = ft_prepare_neighbours(cfg_neighb);
                cfg.neighbours = neighbs;

                % STATS  
                stat = ft_timelockstatistics(cfg, B_all{:}, B_contr{:}); %R_all

%                 cd(alpha_folder); 
                cd(prob_folder);
                save('stat2_regr_ER1_all_semiLog_nr_22_sn','stat','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_regr_ER2semiLog_hp_22','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s
               
                % Check if there are significant clusters after correction
                stat.posclusters(1); stat.negclusters(1);
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END