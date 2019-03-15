%% Clean electrophysiological data from behavioral bad trials
%=========================================================================%
% AUTHOR: Bianca Trovo
% DATE: October/November 2018
% EXPERIMENT: Timelimit_2018

%{
    
    SCOPE: Create indexes for the trials in which participants
    lifted their finger too early or too late. The remaining trials are the
    one after bad artifact rejection (done with
    *myft_preProcess_and_clean_data* (Aaron, 2011).

    OUTPUT: GOOD_BEHAV newgood_resps outliers good_resps_cond  goodrespsall idx_goodtrls idx_allbadtrls whichEarly whichLate conds
    
    Important: you will specially need the variable *idx_allbadtrls*

    WARNING: *isoutlier* matlab function works only on version 2018 onward


%}

% FIX ME: outliers part not implemented or not working yet. Not sure if we
% need to divide by condition since we do this in the averaging part. 
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

%% Load preprocessed files

if subjnum== 1 || subjnum== 18 || subjnum== 19 || subjnum== 20 || subjnum== 21 || subjnum== 22
    load(sprintf('TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp',subjnum))
else
    load(sprintf('TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp',subjnum))
end

%% Create indexes for all trials together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% conds=[];resps=[];clockstarts=[];timediff=[];waitingtimes=[];resptimes=[];

condinf = [TRIALS.cond]';
condinf(condinf==32) = Inf; %replace all the 32 with Inf
conds= condinf; %to not get confused
uconds= unique(conds');
ntrials= length([TRIALS.rt]');
    
    % From samples to seconds
    resps= [TRIALS.rt]';
    clockstarts= [TRIALS.t0]';
    timediff= resps - clockstarts;
    waitingtimes= timediff/500; %divided by sampling rate (500Hz)
    resptimes= waitingtimes-3.0; % from Zafer's code, we know for sure it's 3 sec exaclty.

    % Main variable that we want
    idx_allbadtrls=find(resptimes <0|resptimes > conds); %<=0.2 replaced with 0

    % Extra variables 
    idx_tooEarly=find(resptimes <0); %<=0.2 replaced with 0
    idx_tooLate= find(resptimes > conds);
    idx_goodtrls= find(resptimes >0 & resptimes < conds);
    goodrespsall= resptimes(idx_goodtrls); % in seconds
    whichEarly= resptimes(idx_tooEarly);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Compute percentage just to have an estimation
        perc_early= numel(whichEarly)*100/numel(resptimes);
        whichLate= resptimes(idx_tooLate);
        perc_late= numel(whichLate)*100/numel(resptimes);
        % Display them ?
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% newconds= conds(idx_goodtrls);         
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create indexes divided by condition and remove OUTLIERS [ still not functioning or not implemented ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trials
% good_resps_cond={}; TF={}; outliers={};newgood_resps={};
%     
%     for condi = 1:length(uconds)
%         
%         good_resps_cond{condi}= goodrespsall(find(newconds == uconds(condi))); % BEFORE removing outliers 
%         TF{condi} = isoutlier(good_resps_cond{condi}, 'median');
%         outliers{condi}= good_resps_cond{condi}(TF{condi});
%         newgood_resps{condi}= good_resps_cond{condi}(~TF{condi}); % WITHOUT OUTLIERS
%         
%     end
% 
% % Indexes
% idx_goodxcond={};idxTF={};idx_outl={};idx_newgoodresps={};
% 
%     for condi = 1:length(uconds)
%         
%         idx_goodxcond{condi}=idx_goodtrls(find(newconds == uconds(condi))); % BEFORE removing outliers
%         idxTF{condi} = idx_goodxcond{condi}(TF{condi});
% %         idx_outl{condi}= idx_goodxcond{condi}(idxTF{condi});
% %         idx_newgoodresps{condi}=idx_goodxcond{condi}(~idxTF{condi});
%         
%     end

% ultimate trials per condition without outliers
% good_behav_rt{condi}= idx_outl(find(newconds == uconds(condi))); 

%% Save stuff in the right folder

% Create the folder if it doesn't exist already.
if input('Save BEHAVIORAL INDEXES results? RISK OF OVERWRITING  (1/0) ... ')
    
    behav_folder= [current_subj_folder, '/Behavioral'];
         if ~exist(fullfile(behav_folder)); mkdir(fullfile(behav_folder)); end;
    cd(behav_folder);

%     save GOOD_BEHAV newgood_resps outliers good_resps_cond  goodrespsall idx_goodtrls idx_allbadtrls whichEarly whichLate conds
end

%%
disp(['Done with subject ' num2str(subjnum)])

%% End of the script