%% Compute average & variability within subject
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on November 2018
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the timelocked average (mean RP/RF) and the EEG
    variability (std RP) prior to self-initiated action (event: motor response). For
    the latter we follow Khalighinejad et al. 2018.

    OUTPUT: avg{condi} in avg , across_stdev{condi}, within_std{condi} in
    std.

    FIXME: loop for within-trial variability (it doesn't compute for all
    channels).
 

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

%% Load file 

load([current_subj_folder,'/Behavioral/GOOD_BEHAV.mat']);

%% Compute average here (*timelockanalysis*)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up indexes for getting only the good trials
good_trls = setdiff([1:length(DATA_REJ_INTERP.trial)],idx_allbadtrls);
isequal(idx_goodtrls',good_trls)
% redundant but we redo it just in case
cond= [TRIALS.cond]; %we put all the conditions in a row
cond(cond==32) = Inf;
un_conds = unique(cond);
newcond= cond(good_trls);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-preprocessing data for further analyses 
    cfg=[];
    cfg.trials = good_trls;
    % cfg.lpfilter= 'yes';
    %cfg.lpfreq = 30; % alternatively filter at 20Hz (2 Hz!!);
    cfg.demean='yes';
    % cfg.baselinewindow = [-0.005 0.005]; % Khalighnejad/Haggard's baseline
    DATA= ft_preprocessing(cfg,DATA_REJ_INTERP);
    % DATA_bl = ft_preprocessing(cfg,DATA_REJ_INTERP); % In case you apply
    % baseline 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        avg=[];
        for condi = 1:length(un_conds)
            
            cfg=[];
            cfg.trials= find(newcond == un_conds(condi));
            avg{condi} = ft_timelockanalysis(cfg,DATA);
        %     avg{condi} = ft_timelockanalysis(cfg,DATA_bl);
        %     mean_avg{condi} = avg{condi}.avg; % ???
        
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

check= isequal(avg{1}.avg,avg{2}.avg,avg{3}.avg,avg{4}.avg,avg{5}.avg)

%% Visualization averages timeseries 
% It can be improved with other plot functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timeseries

figure
cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
% cfg.preproc.lpfilter= 'yes';
% cfg.preproc.lpfreq= 35;
cfg.linewidth = 2;
ft_multiplotER(cfg,avg{1},avg{2},avg{3},avg{4},avg{5})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topography? 

%% Save averages timeseries

% Create the folder if it doesn't exist already.
if input('Save TIMELOCKED AVERAGES results? RISK OF OVERWRITING  (1/0) ... ')
    
    timeseries_folder= [current_subj_folder, '/Timeseries'];
         if ~exist(fullfile(timeseries_folder)); mkdir(fullfile(timeseries_folder)); end;
    cd(timeseries_folder);

    save RP_avg avg % choose a proper name and do not change it
    save RP_avg_bl avg % if baseline corrected
    
end

%% Compute variability here (cfr. Khalighinejad et al. 2018)
% You need to use the output of average computation above (=from
% timelockanalysis). Because the fieldtrip function doesn't read in the
% 'var' alone, we need to use the 'avg' denomination, so we take the variability from
% the overall file that contains also the average information (avg.avg
% would be the mean only).

%% Across-trial variability
% WARNING: needs to be baseline corrected to replicate the same findings
across_stdev=[];

    for condi = 1:length(un_conds)
        
        across_stdev{condi} = avg{condi}; % we rename it to avoid confusion
        across_stdev{condi}.avg = sqrt(across_stdev{condi}.var); % we call it avg but it's std
        
    end
 
%% Visualization variability 

hold on;
figure
cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 2;
ft_multiplotER(cfg,across_stdev{1},across_stdev{2},across_stdev{3},across_stdev{4},across_stdev{5})

%% Within-trial variability (check why it wasn't working for multiple channels)

within_std = [];

    for condi = 1:length(un_conds)
        
        good_trls= find(newcond == un_conds(condi));
        for i=1:length(good_trls) % length of goodtrls
            
            within_std{condi}(i,:) = movstd(DATA_bl.trial{i}(20,:),25); %50
            
        end
    end


%% Visualization variability 
% Change the time axis to put seconds instead of samples
tAx = [0:2000]./SR - 3;

for condi = 1:length(un_conds)
    f3(condi)= figure(condi)
    plot(mean(within_std{condi}))
    hold on
end

%% Save variability timeseries

% Create the folder if it doesn't exist already.
if input('Save TIMELOCKED VARIABILITY results? RISK OF OVERWRITING  (1/0) ... ')
    
    timeseries_folder= [current_subj_folder, '/Timeseries'];
         if ~exist(fullfile(timeseries_folder)); mkdir(fullfile(timeseries_folder)); end;
    cd(timeseries_folder);

    save RP_std across_stdev within_std % choose a proper name and do not change it
    
end

% save std_EEG across_stdev within_std 

%% %% End (for now)
disp(['END of the script for subj ' int2str(subjnum)])
clear all; close all;
