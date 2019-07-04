%% Compute time-series within subject
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: re_created on May 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the movement-timelocked ERP average within subject and
    between subject (grandaverage)for all channels for all time points
    within the window -3s 0s.

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
    
%% Load preprocessed files, clean from behavioural artifacts, postprocess

for subi=1:nSubjs; %nSubjs
    
    cd([data_Path, sprintf('/subj%02d', subi)])
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-preprocessing data for further analyses
    cfg=[];
    cfg.trials = good_trls;
    cfg.lpfilter= 'yes';
    cfg.lpfreq = 40; % alternatively filter at 20Hz (2 Hz!!);
    cfg.demean='yes';
    % cfg.baselinewindow = [-0.005 0.005]; % Khalighnejad/Haggard's baseline
%     DATA_CLEAN= ft_selectdata(cfg,DATA_REJ_INTERP);
    DATA_CLEAN= ft_preprocessing(cfg,DATA_REJ_INTERP);
    % DATA_bl = ft_preprocessing(cfg,DATA_REJ_INTERP); % If you apply baseline
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % avg_cond
    
    avg_cond=[];
    
    for condi = 1:length(un_conds)
        
        cfg=[];
        cfg.trials= find(newcond == un_conds(condi));
        avg_cond{condi} = ft_timelockanalysis(cfg,DATA_CLEAN);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % avg_trl
    
    avg_trl=[];
   
    cfg=[];
    cfg.keeptrials  = 'yes';
    avg_trl = ft_timelockanalysis(cfg,DATA_CLEAN);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % avg_one : you will need it for the template for the cluster test
    
    avg_one=[];
    
    cfg=[];
    avg_one = ft_timelockanalysis(cfg,DATA_CLEAN);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % avg_condTrl : you will need it for the fancy plot
    
    avg_condTrl=[];
    
    for condi = 1:length(un_conds)
        
        cfg=[];
        cfg.trials= find(newcond == un_conds(condi));
        cfg.keeptrials  = 'yes';
        avg_condTrl{condi} = ft_timelockanalysis(cfg,DATA_CLEAN);
        
    end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        
    %% Save averages time-series
    
    % Create the folder if it doesn't exist already
    
    timeseries_folder= [results_Path, '/Timeseries']; % it can be also current_subj_folder
    if ~exist(fullfile(timeseries_folder)); mkdir(fullfile(timeseries_folder)); end;
    cd(timeseries_folder);
    
    % avg_cond
    filename= [sprintf('subj%02d_TimeS_bytrial', subi)]; % add one if all trials mixed by condition
    save(filename,'avg_cond','-v7.3');
    % avg_trl
    filename= [sprintf('subj%02d_TimeS_bytrial', subi)]; % add one if all trials mixed by condition
    save(filename,'avg_trl','-v7.3');
    % avg_one
    filename= [sprintf('subj%02d_TimeS_bytrial', subi)]; % add one if all trials mixed by condition
    save(filename,'avg_one','-v7.3');
    % avg_condTrl
    filename= [sprintf('subj%02d_TimeS_bytrial', subi)]; % add one if all trials mixed by condition
    save(filename,'avg_condTrl','-v7.3');
    
    disp(['Subject ' num2str(subi) ' done']);

end % of the whole loop across nSubjs

% filename= [sprintf('subj%02d_usefulinfo', subi)]; % add one if all trials mixed by condition
% save(filename,'newcond','-v7.3');
% disp(['Subject ' num2str(subi) ' done']);
% 
% end

cd([results_Path, '/Timeseries']);
%disp('END of TIME-SERIES ANALYSIS by condition per subject');
disp('END of TIME-SERIES ANALYSIS per subject');

%% GRANDAVERAGE 
nSubjs= 22;
%~~~ Other
OKsubjs= [3 6 7 8 10 13 15 17 18 19 20 21]; % done with the 1µV method assessment 

cd(pwd);
for subi=1:nSubjs;
    
    fname_ER= sprintf('subj%02d_TimeS_one',subi);
    pickupER(subi) = load(fname_ER);
   
end
% 
% save pickupER_cond pickupER;

%%
TimeSmatrix=[];

for subi= 1:nSubjs; 
 
    for k= 1:5
        TimeSmatrix{subi,k}= pickupER(subi).avg_cond{k}; % or avg_condTrl or avg_trl or avg_one
    end

end

save TimeSmatrix TimeSmatrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Grand_ER=[]; 

for k= 1:5
    
    cfg=[];
%     cfg.keepindividual= 'yes';
%     cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
    Grand_ER{k}= ft_timelockgrandaverage(cfg,TimeSmatrix{:,k});

end

save Grand_ER_Ind Grand_ER

%% END of the script

fprintf('\n END: DATA SAVED \n')