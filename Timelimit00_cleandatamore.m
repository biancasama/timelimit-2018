%% 

% FIXME: Add option for chosing between EEG and MEG
% the script is getting too long. do we really need to compute all the
% variables for 
%% START: initialize path 21st Nov

clearvars;
% Define some paths 
if strcmp(computer, 'MACI64')% on my laptop
    script_Path= '/Volumes/USB_DISK/TIMELIMIT_backup/SCRIPTS_ANALYSES/MEEG'; % here you find all the scripts for preprocessing/analysing MEEG data
    data_Path = '/Volumes/USB_DISK/TIMELIMIT_backup/MEEG_fif_files'; % = parent_folder: all the raw data are stored here.
    ft_Path = '/Users/bt_neurospin/matlab/FIELDTRIP/fieldtrip-20170405'; % Fieldtrip tools 
    tool_Path= '/Users/bt_neurospin/matlab'; % other useful functions for matlab (written by Aaron)
end

% add some paths
addpath(genpath(script_Path)); % general pre-proc path
addpath(genpath(tool_Path));
addpath(ft_Path); ft_defaults; % start up fieldtrip [NEW]
addpath([ft_Path '/engine']); % start up FT engines [NEW] NOTE: what is this actually doing? 

%% Choose the subject
prompt={'Which participant do you want to look at?'};
name='Subject number';
numlines=1;
answer=inputdlg(prompt,name,numlines);
subjnum= str2double(answer);

%% Move to the right folder/path 
% 
parent_folder='/Volumes/USB_DISK/TIMELIMIT_backup/MEEG_fif_files'; 
% parent_folder= '/Volumes/LaCie/128_usb_BACKUP/Project_Timelimit/fif_files_timelimit';
subj_folders = dir(fullfile(parent_folder, 'subj*'));
current_subj_folder = fullfile(parent_folder, subj_folders(subjnum).name);
cd(current_subj_folder);

%% Load data preprocessed (SHOWS PROBLEM ON MATLAB 2018) 

% if subjnum== 1
%     load('TimeLimit_v2_Resp_subj01_EEG_clean_concat_rej_interp.mat');
% else
    filename = sprintf('Time_Limit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp.mat',subjnum); %'subj%02d/avg_EEG.mat'
    load(filename);
% end

%%

conds=[];resps=[];clockstarts=[];timediff=[];waitingtimes=[];resptimes=[];
% AlltooEarly = []; AlltooLate = [];


    condinf = [TRIALS.cond]';
    condinf(condinf==32) = Inf; %replace all the 32 with Inf (Elie)
    conds= condinf; %to not get confused
    uconds= unique(conds');
    ntrials= length([TRIALS.rt]');
    
    %Added 16th November
    resps= [TRIALS.rt]';
    clockstarts= [TRIALS.t0]';
    timediff= [resps - clockstarts];
    waitingtimes= timediff/500; %divided by sampling rate (500Hz)
    resptimes= waitingtimes-3.0; %DO WE KNOW IF IT'S ALWAYS 3 SEC EXACTLY?
    idx_tooEarly=find(resptimes <0); %<=0.2 replaced with 0
    idx_tooLate= find(resptimes > conds);
    idx_allbadtrls=find(resptimes <0|resptimes > conds); %<=0.2 replaced with 0
    idx_goodtrls= find(resptimes >0 & resptimes < conds);
    goodrespsall= resptimes(idx_goodtrls)
    whichEarly= resptimes(idx_tooEarly)
    perc_early= numel(whichEarly)*100/numel(resptimes)
%     disp([num2str(perc_early)' of the cleaned trials were too early responses'])
    whichLate= resptimes(idx_tooLate)
    perc_late= numel(whichLate)*100/numel(resptimes)
    
    newconds= conds(idx_goodtrls); 
   
    good_resps_cond={}; TF={}; outliers={};newgood_resps={};
    
    for condi = 1:length(uconds)
        good_resps_cond{condi}= goodrespsall(find(newconds == uconds(condi)))
%         idx_goodxcond{condi}=idx_goodtrls(find(newconds == uconds(condi)))
        TF{condi} = isoutlier(good_resps_cond{condi}, 'median');
%         idxTF{condi} = idx_goodxcond{condi}(TF{condi})
        outliers{condi}= good_resps_cond{condi}(TF{condi});
%         idx_outl{condi}= idx_goodxcond{condi}(idxTF{condi})
        newgood_resps{condi}= good_resps_cond{condi}(~TF{condi})
%         idx_newgoodresps{condi}=idx_goodxcond{condi}(~idxTF{condi})
    end
    
    idx_goodxcond={};idxTF={};idx_outl={};idx_newgoodresps={};
    for condi = 1:length(uconds)
        
        idx_goodxcond{condi}=idx_goodtrls(find(newconds == uconds(condi)));
        
        idxTF{condi} = idx_goodxcond{condi}(TF{condi});
       
%         idx_outl{condi}= idx_goodxcond{condi}(idxTF{condi})
%         
%         idx_newgoodresps{condi}=idx_goodxcond{condi}(~idxTF{condi})
        
    end
    
% % 
%     TF = isoutlier(A, METHOD)
%     [TF, LTHRESH, UTHRESH, CENTER] = isoutlier(...) 
% TF = isoutlier(good_resps_cond{3}, 'median')'
% outliers= good_resps_cond{3}(TF)
% newgood_resps=good_resps_all(~TF)

% good_behav_rt{condi}= idx_outl(find(newconds == uconds(condi)))

%% Save stuff in the right folder

% Create the folder if it doesn't exist already.
status = mkdir(current_subj_folder,'Behavioral');
behav_folder= fullfile(current_subj_folder,'Behavioral');
cd(behav_folder);

% save GOOD_BEHAV newgood_resps good_resps_cond goodrespsall outliers whichEarly whichLate conds
save GOOD_BEHAV newgood_resps outliers good_resps_cond  goodrespsall idx_goodtrls idx_allbadtrls whichEarly whichLate conds

%%
disp(['Done with subject ' num2str(subjnum)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to by condition
% 
%     idx_cond1=[];idx_cond2=[];idx_cond3=[];idx_cond4=[];idx_cond5=[];
%     for condi = 1:length(uconds)
%         
%         idx_cond1= find([TRIALS.rt] & [TRIALS.cond] == uconds(1));
%         idx_cond2= find([TRIALS.rt] & [TRIALS.cond] == uconds(2));
%         idx_cond3= find([TRIALS.rt] & [TRIALS.cond] == uconds(3));
%         idx_cond4= find([TRIALS.rt] & [TRIALS.cond] == uconds(4));
%         idx_cond5= find([TRIALS.rt] & [TRIALS.cond] == uconds(5));
%     end
% 
% %     idx_allconds= [idx_cond1; idx_cond2; idx_cond3; idx_cond4; idx_cond5]
%     % SORT EVERYTHING IN 5 CONDITIONS 
%     
%     resps= [];
%     resps= [TRIALS.rt]';
%     resp_cond1= []; resp_cond2= []; resp_cond3= []; resp_cond4= []; resp_cond5= [];
%         resp_cond1= resps(idx_cond1);
%         resp_cond2= resps(idx_cond2);
%         resp_cond3= resps(idx_cond3);
%         resp_cond4= resps(idx_cond4);
%         resp_cond5= resps(idx_cond5);
%     
%     RESPS{1}= resp_cond1; RESPS{2}= resp_cond2; RESPS{3}= resp_cond3;
%     RESPS{4}= resp_cond4; RESPS{5}= resp_cond5;
%     
%     clockstarts= [];
%     clockstarts= [TRIALS.t0]';
%     start_cond1= []; start_cond2= []; start_cond3= []; start_cond4= []; start_cond5= [];
%         start_cond1= clockstarts(idx_cond1);
%         start_cond2= clockstarts(idx_cond2);
%         start_cond3= clockstarts(idx_cond3);
%         start_cond4= clockstarts(idx_cond4);
%         start_cond5= clockstarts(idx_cond5);
%         
%     START{1}= start_cond1; START{2}= start_cond2; START{3}= start_cond3; 
%     START{4}= start_cond4; START{5}= start_cond5;
% 
% %     timediff= [resps - clockstarts];
% TIMEDIFF=[]
%     for condi = 1:length(uconds)
%         TIMEDIFF{condi}= [RESPS{condi}-START{condi}];
%     end
% %     
% %      waitingtimes= timediff/500; %divided by sampling rate (500Hz)
%   WAITINGTIME= []  
%   for condi = 1:length(uconds)
%       WAITINGTIME{condi}= TIMEDIFF{condi}/500
%   end
%                     
% %     resptimes= waitingtimes-3.0; %DO WE KNOW IF IT'S ALWAYS 3 SEC EXACTLY?
%     REALRESPTIME=[]
%     for condi = 1:length(uconds)
%         REALRESPTIME{condi}= WAITINGTIME{condi}-3.0
%     end
%                         
% %     tooEarly=find(resptimes <0); %<=0.2 replaced with 0
% %     tooLate= find(resptimes > conds);
% %     allbad_trls=find(resptimes <0|resptimes > conds); %<=0.2 replaced with 0
% %     good_trls= find(resptimes >0 & resptimes < conds);
%     
%     idxEARLY=[];idxLATE=[];idxBAD=[];idxGOOD=[];
%     for condi = 1:length(uconds)
%         idxEARLY{condi}=find(REALRESPTIME{condi} <0); %<=0.2 replaced with 0
%         idxLATE{condi}= find(REALRESPTIME{condi} > uconds(condi));
%         idxBAD{condi}=find(REALRESPTIME{condi} <0|REALRESPTIME{condi} > uconds(condi)); %<=0.2 replaced with 0
%         idxGOOD{condi}= find(REALRESPTIME{condi} >0 & REALRESPTIME{condi} < uconds(condi));
%     end
%     
%     % no need anymore
% % whichBAD= conds(allbad_trls); % which conditions do the bad trials come from?
% % whichGOOD= conds(good_trls); % added 5th Nov.
% % EarlyRT_trigg= resptimes(tooEarly);
% % LateRT_trigg= resptimes(tooLate);
% 
% 
% 
% % crucial part PROBLEMATIC AFTER THE CELL THING
% finalgood=resptimes(good_trls);
% resptimes(idx_allbadtrls)=[];
% 
% % for condi = 1:length(uconds)
% % TOKEEP= REALRESPTIME{condi}(idxGOOD{condi})
% 
% 
% % if statement: if yes, save variable
% isequal(resptimes,finalgood)
% % add disp YES if 1 NO if 0
% Clean_TRLS= resptimes;
% disp([num2str(numel(idx_allbadtrls)), ' trials have been removed from subject ' num2str(subjnum) '; ' num2str(numel(Clean_TRLS)) ' trials remaining']); % just to confirm
% 
% %% save stuff
% 
% save BADcleanedtrls idx_allbadtrls whichBAD 
% save Clean_eegTRLS Clean_TRLS good_trls whichGOOD
% 
% %% end of the script 