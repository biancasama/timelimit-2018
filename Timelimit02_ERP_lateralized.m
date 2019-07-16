%% Compute time-series within subject
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: July 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: 

    OUTPUT: avg_cond, avg_trl, avg_one, avg_condTrl;  TimeSmatrix,
    Grand_ER.

    FIXME: **

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

clear LevelAnalysis name numlines prompt subj_folders; 

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
    
%% Load ***

load TimeSmatrix;
load Grand_ER_Ind;

%% Compute LRP

load TimeSmatrix; % all conditions in one, trials averaged

for subi= 1:nSubjs;
    for k= 1:5
        cfg=[];
        cfg.channel= 'EEG028';
        TimeS_C3{subi,k}= ft_selectdata(cfg,TimeSmatrix{subi,k});
    end
end

for subi= 1:nSubjs;
    for k= 1:5
        cfg=[];
        cfg.channel= 'EEG032';
        TimeS_C4{subi,k}= ft_selectdata(cfg,TimeSmatrix{subi,k});
    end
end

templateLRP= TimeS_C3; % fake Fieldtrip structure for inserting correlation values
save templateLRP templateLRP;

LRPmatrix={};

for subi= 1:nSubjs; 
 
    for k= 1:5
        %templateER{subi,k}.avg= pickupER(subi).avg_cond{k}.avg(28,:)-pickupER(subi).avg_cond{k}.avg(32,:); % or avg_condTrl or avg_trl or avg_one
        templateLRP{subi,k}.avg= TimeS_C3{subi,k}.avg-TimeS_C4{subi,k}.avg; % or avg_condTrl or avg_trl or avg_one
        LRPmatrix= templateLRP;
    end

end

save LRPmatrix LRPmatrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % for subi= 1:nSubjs;
% for k= 1:5
%     cfg=[];
%     cfg.channel= 'EEG028';
%     Grand_C3{k}= ft_selectdata(cfg,Grand_ER{k});
% end
% % end
% for k= 1:5
%     cfg=[];
%     cfg.channel= 'EEG032';
%     Grand_C4{k}= ft_selectdata(cfg,Grand_ER{k});
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% GrandLRP= [];
% % Then take the difference of the averages using ft_math
% for k= 1:5
%     cfg  = [];
%     cfg.operation = 'subtract';
%     cfg.parameter = 'avg';
%     % cfg.matrix= GrandLRP; 
%     GrandLRP{k}= ft_math(cfg, Grand_C3{k}, Grand_C4{k});
% end

GrandLRP={}; 

for k= 1:5
    
    cfg=[];
    cfg.keepindividual= 'yes';
%     cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
    GrandLRP{k}= ft_timelockgrandaverage(cfg,LRPmatrix{:,k});

end

save GrandLRP GrandLRP;
save GrandLRP_Ind GrandLRP;

%% by trial
pickupLRP=[];
for subi= 1:nSubjs
    for i= 1: length(pickupER(subi).avg_trl.trial(:,1,1))
        pickupLRP(subi).avg_trl.trial=pickupER(subi).avg_trl.trial(:,28,:)- pickupER(subi).avg_trl.trial(:,32,:);
    end
end

save pickupLRP pickupLRP -v7.3;



%% END of the script

fprintf('\n END: DATA SAVED \n')