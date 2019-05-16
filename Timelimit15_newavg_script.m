%% Compute time-series within subject
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: re_created on May 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the 

    OUTPUT: avg{condi} in 

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

%% Load preprocessed files or 

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
    % DATA_bl = ft_preprocessing(cfg,DATA_REJ_INTERP); % In case you apply
    % baseline
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    avg_trl=[];
    %         for condi = 1:length(un_conds)
    
    cfg=[];
    %             cfg.trials= find(newcond == un_conds(condi));
    cfg.keeptrials  = 'yes';
    avg_trl = ft_timelockanalysis(cfg,DATA_CLEAN);
    %             avg{condi} = ft_timelockanalysis(cfg,DATA_CLEAN);
    
    %         end
   
        
    % Save averages time-series
    
    % Create the folder if it doesn't exist already
    
    timeseries_folder= [results_Path, '/Timeseries']; % it can be also current_subj_folder
    if ~exist(fullfile(timeseries_folder)); mkdir(fullfile(timeseries_folder)); end;
    cd(timeseries_folder);
    
    filename= [sprintf('subj%02d_TimeS_bytrial', subi)]; % add one if all trials mixed by condition
    save(filename,'avg_trl','-v7.3');
    
    disp(['Subject ' num2str(subi) ' done']);

end

% filename= [sprintf('subj%02d_usefulinfo', subi)]; % add one if all trials mixed by condition
% save(filename,'newcond','-v7.3');
% disp(['Subject ' num2str(subi) ' done']);
% 
% end

cd([results_Path, '/Timeseries']);
%disp('END of TIME-SERIES ANALYSIS by condition per subject');
disp('END of TIME-SERIES ANALYSIS per subject');

%% Grandaverage
nSubjs= 22;
%~~~ Other
OKsubjs= [3 6 7 8 10 13 15 17 18 19 20 21];

cd(pwd);
for subi=1:nSubjs;
    
    fname_TimeS= sprintf('subj%02d_TimeS_one',subi);
    pickupTimeS(subi) = load(fname_TimeS);
   
end


% %
% TimeSmatrix=[];
% for subi= 1:nSubjs; %nGoodSubjects
%     
%     for k= 1:5
%         TimeSmatrix{subi,k}= pickupTimeS(subi).TFR{k}; %pickupSub(i).avg_EEG{k}
%     end
% end
% 
% save TimeSmatrix TimeSmatrix
% 
% Grand_TFR=[]; 
% for k= 1:5
%     
% %     cfg=[];
% %   
% %     cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
%     Grand_TFR{k}= ft_freqgrandaverage([],TimeSmatrix{:,k});
% %     Grand_Freq_ROI{k}= mean(Grand_Freq{k}.powspctrm,1);
%     
% end
% 
% save Grand_TFR Grand_TFR
% 
% cfg=[];
% cfg.baseline= [-2 -1.5];
% cfg.baselinetype= 'relative';
% % cfg.zlim= [-3e-25 3e-25];
% cfg.xlim= [-1 -0];
% % cfg.ylim= [8 30];
% % cfg.ylim= [15 30];
% cfg.layout = 'eeg_64_NM20884N.lay';
% cfg.linewidth = 2;
% cfg.showlabels= 'yes';
% figure
% ft_multiplotTFR(cfg,Grand_TFR);%Grand_TFR{2},Grand_TFR{3},Grand_TFR{4},Grand_TFR{5}
% % ft_topoplotTFR(cfg,Grand_TFR{2});