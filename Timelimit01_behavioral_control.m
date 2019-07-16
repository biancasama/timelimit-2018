%% Behavioral analysis based on MEEG data
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on July 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the 

    OUTPUT: 

    FIXME: SEM for medians

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
    
regression_folder= [results_Path, '/Regressions']; % it can be also current_subj_folder
    if ~exist(fullfile(regression_folder)); mkdir(fullfile(regression_folder)); end;
    
%% Load preprocessed files and exclude behavioural artifacts 

for subi=1:nSubjs; %nSubjs
    
        cd([data_Path, sprintf('/subj%02d', subi)]);
%     current_subj_folder= fullfile(data_Path, subj_folders(subi).name);
%     cd(current_subj_folder);
    
    if subi== 1 || subi== 18 || subi== 19 || subi== 20 || subi== 21 || subi== 22
        load(sprintf('TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp',subi))
    else
        load(sprintf('TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp',subi))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setting up indexes for getting only the good trials
    
    [idx_goodxcond,idx_goodtrls, idx_allbadtrls]= BTmy_cleandatamore(TRIALS);
    
    good_trls = setdiff([1:length(DATA_REJ_INTERP.trial)],idx_allbadtrls);
    
    if isequal(idx_goodtrls',good_trls)==1; disp('YES'); else disp('NO'); end;
    
    % redundant but we redo it just in case
    cond= [TRIALS.cond]; %we put all the conditions in a row
    cond(cond==32) = Inf;
    un_conds = unique(cond);
    newcond= cond(good_trls);
    % add control
    idx_cond2keep = find(newcond < Inf); % either > 2 or < Inf
    contrl_cond= newcond(idx_cond2keep);
    un_conds2 = unique(contrl_cond);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-preprocessing data for further analyses
    
    cfg=[];
    cfg.trials = idx_cond2keep;
    
    DATA_CLEAN= ft_selectdata(cfg,DATA_REJ_INTERP);
    
    resps= [TRIALS.rt];
    clockstarts= [TRIALS.t0];
    TRIALS_CLEAN= resps(idx_cond2keep);
    CLOCKST= clockstarts(idx_cond2keep);
    
    if isequal(length(DATA_CLEAN.trial), length(TRIALS_CLEAN), length(CLOCKST)); disp(['CORRECT: M/EEG data of subj ' num2str(subi) ' matches size BEHAV data']); else disp(['INCORRECT:M/EEG data of subj ' num2str(subi) 'does not matche size BEHAV data']); end;
    
  
    TIMEDIFF=  [TRIALS_CLEAN- CLOCKST];
    WAITTIMES= TIMEDIFF/500; %divided by sampling rate (500Hz)
    RESPTIMES=  WAITTIMES-3.0; % from Zafer's code, we know for sure it's 3 sec exaclty.
    
    % Sorted
    contr_resps_cond={}; 
    
    for condi = 1:length(un_conds2)
        
        contr_resps_cond{condi}= RESPTIMES(find(contrl_cond == un_conds2(condi))); % BEFORE removing outliers
        
    end
   
    behavioral_folder= [results_Path, '/Behaviour']; % it can be also current_subj_folder
    if ~exist(fullfile(behavioral_folder)); mkdir(fullfile(behavioral_folder)); end;
    cd(behavioral_folder);
    
    filename= [sprintf('subj%02d_WaitingTimesContr', subi)]; % add one if all trials mixed by condition
    save(filename,'RESPTIMES','contr_resps_cond','contrl_cond','idx_cond2keep');
    
    disp(['End of subj ' num2str(subi)]);
    
end

%% Load/group behavioral data, compute log and normalize resp times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
behavioral_folder2= [behavioral_folder, '/BehaviourContr']; % it can be also current_subj_folder
    if ~exist(fullfile(behavioral_folder2)); mkdir(fullfile(behavioral_folder2)); end;
cd(behavioral_folder2);
    
% load in an unique structure 
for subi=1:nSubjs
    %
    %         cd(behavioral_folder2);
    fname_BehavData= sprintf('subj%02d_WaitingTimesContr',subi);
    pickupBehav(subi) = load(fname_BehavData);
    
end

% load in an unique structure 
for subi=1:nSubjs
    %
    %         cd(behavioral_folder2);
    fname_Cond= sprintf('subj%02d_WaitingTimesContr',subi);
    pickupCond(subi) = load(fname_Cond);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for subi=1:nSubjs;
    pickupBehav(subi).LogRESPS = log(pickupBehav(subi).RESPTIMES);
    
    for condi= 1:length(un_conds2)
        pickupBehav(subi).LogByConds{condi} = log(pickupBehav(subi).contr_resps_cond{condi});
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalization with method 'z-score'
for subi=1:nSubjs;
    
    pickupBehav(subi).normRESPS = normalize(pickupBehav(subi).RESPTIMES);
    
    for condi= 1:length(un_conds2)
        pickupBehav(subi).normRESPCond{condi} = normalize(pickupBehav(subi).contr_resps_cond{condi});
    end
    
end

save pickupBehavContr pickupBehav;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% just to see
ALL_BEHAV= [pickupBehav(subi).RESPTIMES'  pickupBehav(subi).LogRESPS'];

%[pickupBehav(subi).good_resps_cond{1}' pickupBehav(subi).good_resps_cond{2}' pickupBehav(subi).good_resps_cond{3}' pickupBehav(subi).good_resps_cond{4}' pickupBehav(subi).good_resps_cond{5}']

%% Descriptive stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
behavStats= struct('mWT',[],'mdWT', [],'stdWT',[], 'semWT',[], 'minWT',[], 'maxWT',[]); 
for subi= 1: nSubjs
    for condi= 1:5
        
    behavStats.mWT(subi,condi)= nanmean(pickupBehav(subi).good_resps_cond{condi}); % mean
    behavStats.mdWT(subi,condi)= nanmedian(pickupBehav(subi).good_resps_cond{condi}); % median
    behavStats.stdWT(subi,condi)= nanstd(pickupBehav(subi).good_resps_cond{condi}); % standard deviation
    behavStats.semWT(subi,condi)= sem(pickupBehav(subi).good_resps_cond{condi}); % Standard Error of the Mean
    behavStats.minWT(subi,condi)= nanmin(pickupBehav(subi).good_resps_cond{condi}); % minimum
    behavStats.maxWT(subi,condi)= nanmax(pickupBehav(subi).good_resps_cond{condi}); % maximum
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interquartile range for median values

for condi= 1:5
%     
    pd(condi) = fitdist(behavStats.mdWT(:,condi),'Normal'); % probability distribution
    r(condi) = iqr(pd(condi));  %interquartile range values (in matlab outputs appear as 'r'
    y(condi,:) = icdf(pd(condi),[0.25,0.75]); % Inverse cumulative distribution function 
    
end

IQR= y; % saved independently in DescriptiveStats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% max avg response

for condi=1:5; maxWaits(condi)= max(behavStats.mdWT(:,condi)); semMaxWaits(condi)= sem(behavStats.mdWT(:,condi)); end;
figure; bar(maxWaits);

% min avg response

for condi=1:5; minWaits(condi)= min(behavStats.mdWT(:,condi)); end;
figure; bar(minWaits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in log scale
LogBehavStats= struct('mWT',[],'mdWT', [],'stdWT',[], 'semWT',[], 'minWT',[], 'maxWT',[]); 
for subi= 1: nSubjs
    for condi= 1:5
        
    LogBehavStats.mWT(subi,condi)= nanmean(pickupBehav(subi).LogByConds{condi});
    LogBehavStats.mdWT(subi,condi)= nanmedian(pickupBehav(subi).LogByConds{condi});
    LogBehavStats.stdWT(subi,condi)= nanstd(pickupBehav(subi).LogByConds{condi});
    LogBehavStats.semWT(subi,condi)= sem(pickupBehav(subi).LogByConds{condi});
    LogBehavStats.minWT(subi,condi)= nanmin(pickupBehav(subi).LogByConds{condi});
    LogBehavStats.maxWT(subi,condi)= nanmax(pickupBehav(subi).LogByConds{condi});
   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized (sanity check because mean=0 and std=1)
normStats= struct('mWT',[],'mdWT', [],'stdWT',[], 'semWT',[], 'minWT',[], 'maxWT',[]); 
for subi= 1: nSubjs
    for condi= 1:5
        
    normStats.mWT(subi,condi)= nanmean(pickupBehav(subi).normRESPCond{condi});
    normStats.mdWT(subi,condi)= nanmedian(pickupBehav(subi).normRESPCond{condi});
    normStats.stdWT(subi,condi)= nanstd(pickupBehav(subi).normRESPCond{condi});
    normStats.semWT(subi,condi)= sem(pickupBehav(subi).normRESPCond{condi});
    normStats.minWT(subi,condi)= nanmin(pickupBehav(subi).normRESPCond{condi});
    normStats.maxWT(subi,condi)= nanmax(pickupBehav(subi).normRESPCond{condi});
   
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grandaverages
GAVGbehav= [];
for condi= 1:5
    
    GAVGbehav.mWT(condi)= nanmean(behavStats.mWT(:,condi),1);
    GAVGbehav.mdWT(condi)= nanmean(behavStats.mdWT(:,condi),1);
    GAVGbehav.stdWT(condi)= nanmean(behavStats.stdWT(:,condi),1);
    GAVGbehav.semWT(condi)= nanmean(behavStats.semWT(:,condi),1);
    
    GAVGbehav.minWT(condi)= nanmean(behavStats.minWT(:,condi),1);
    GAVGbehav.maxWT(condi)= nanmean(behavStats.maxWT(:,condi),1);
    GAVGbehav.semMaxWT(condi)= sem(behavStats.maxWT(:,condi),1);
    GAVGbehav.semMinWT(condi)= sem(behavStats.minWT(:,condi),1);
    
end

% Log scale

GAVGLogbehav= [];
for condi= 1:5
    
    GAVGLogbehav.mWT(condi)= nanmean(LogBehavStats.mWT(:,condi),1);
    GAVGLogbehav.mdWT(condi)= nanmean(LogBehavStats.mdWT(:,condi),1);
    GAVGLogbehav.stdWT(condi)= nanmean(LogBehavStats.stdWT(:,condi),1);
    GAVGLogbehav.semWT(condi)= nanmean(LogBehavStats.semWT(:,condi),1);
    
    GAVGLogbehav.minWT(condi)= nanmean(LogBehavStats.minWT(:,condi),1);
    GAVGLogbehav.maxWT(condi)= nanmean(LogBehavStats.maxWT(:,condi),1);
    
end

% Save variables 

save DescriptiveStatsContr behavStats  LogBehavStats GAVGbehav GAVGLogbehav IQR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END of the script
