%%

% 
% % for subi= 1:nSubj
% RESPS= []; resps= [];
%     for i= 1:10
%         
%         resps= [TRIALS_CELL{i}.rt]';
%         clockstarts= [TRIALS_CELL{i}.t0]';
%         CLOCKST= [CLOCKST
%             clockstarts];
%         RESPS= [RESPS
%             resps];
%         
%     end
%     save it in each subj folder
%     end
%  
%   
%     TIMEDIFF= RESPS - CLOCKST;
%     waitingtimes= timediff/500; %divided by sampling rate (500Hz)
%     resptimes= waitingtimes-3.0; % from Zafer's code, we know for sure it's 3 sec exaclty.
% 
%     % Main variable that we want
%     idx_allbadtrls=find(resptimes <= 0|resptimes >= conds); %<=0.2 replaced with 0
% 
%     % Extra variables 
%     idx_tooEarly=find(resptimes <= 0); %<=0.2 replaced with 0
%     idx_tooLate= find(resptimes >= conds);
%     idx_goodtrls= find(resptimes > 0 & resptimes < conds);
%     goodrespsall= resptimes(idx_goodtrls); % in seconds
%     whichEarly= resptimes(idx_tooEarly);   

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

for subjnum=1:nSubjs; %nSubjs
    
        cd([data_Path, sprintf('/subj%02d', subjnum)]);
%     current_subj_folder= fullfile(data_Path, subj_folders(subi).name);
%     cd(current_subj_folder);
    
    if subjnum== 1 || subjnum== 18 || subjnum== 19 || subjnum== 20 || subjnum== 21 || subjnum== 22
        load(sprintf('TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp',subjnum))
    else
        load(sprintf('TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp',subjnum))
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
    cfg.trials = find(newcond == un_conds(condi));
    %             cfg.latency = [-1 -.2];
    DATA_CLEAN= ft_selectdata(cfg,DATA_REJ_INTERP);
    
    resps= [TRIALS.rt];
    clockstarts= [TRIALS.t0];
    TRIALS_CLEAN= resps(good_trls);
    CLOCKST= clockstarts(good_trls);
    
    if isequal(length(DATA_CLEAN.trial), length(TRIALS_CLEAN), length(CLOCKST)); disp(['CORRECT: M/EEG data of subj ' num2str(subjnum) ' matches size BEHAV data']); else disp(['INCORRECT:M/EEG data of subj ' num2str(subjnum) 'does not matche size BEHAV data']); end;
    
  
    TIMEDIFF=  [TRIALS_CLEAN- CLOCKST];
    WAITTIMES= TIMEDIFF/500; %divided by sampling rate (500Hz)
    RESPTIMES=  WAITTIMES-3.0; % from Zafer's code, we know for sure it's 3 sec exaclty.
    
    % Sorted
    good_resps_cond={}; 
    
    for condi = 1:length(uconds)
        
        good_resps_cond{condi}= RESPTIMES(find(newcond == un_conds(condi))); % BEFORE removing outliers
        
    end
   
    behavioral_folder= [results_Path, '/Behaviour']; % it can be also current_subj_folder
    if ~exist(fullfile(behavioral_folder)); mkdir(fullfile(behavioral_folder)); end;
    cd(behavioral_folder);
    
    filename= [sprintf('subj%02d_WaitingTimes', subjnum)]; % add one if all trials mixed by condition
    save(filename,'RESPTIMES','good_resps_cond');
    
    disp(['End of subj ' num2str(subjnum)]);
    
end
