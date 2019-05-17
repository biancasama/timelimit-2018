%% Behavioral analysis based on MEEG data
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on April/May 2019
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
%% Load/group behavioral data, compute log and normalize resp times

% load in an unique structure 
for subi=1:nSubjs
    %
    %         cd(behavioral_folder2);
    fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
    pickupBehav(subi) = load(fname_BehavData);
    
end
   
% log tranformation
for subi=1:nSubjs;
    pickupBehav(subi).LogRESPS = log(pickupBehav(subi).RESPTIMES);
    
    for condi= 1:5
        pickupBehav(subi).LogByConds{condi} = log(pickupBehav(subi).good_resps_cond{condi});
    end
    
end

% normalization with method 'z-score'
for subi=1:nSubjs;
    
    pickupBehav(subi).normRESPS = normalize(pickupBehav(subi).RESPTIMES);
    
    for condi= 1:5
        pickupBehav(subi).normRESPCond{condi} = normalize(pickupBehav(subi).good_resps_cond{condi});
    end
    
end

save pickupBehav pickupBehav;

% just to see
ALL_BEHAV= [pickupBehav(subi).RESPTIMES'  pickupBehav(subi).LogRESPS'];

%[pickupBehav(subi).good_resps_cond{1}' pickupBehav(subi).good_resps_cond{2}' pickupBehav(subi).good_resps_cond{3}' pickupBehav(subi).good_resps_cond{4}' pickupBehav(subi).good_resps_cond{5}']

%% descriptive stats

 for subi=1:nSubjs;

%         cd(powerspectra_folder);
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);
        

 end


behavStats= struct('mWT',[],'mdWT', [],'stdWT',[], 'semWT',[], 'minWT',[], 'maxWT',[]); 
for subi= 1: nSubjs
    for condi= 1:5
        
    behavStats.mWT(subi,condi)= nanmean(pickupBehav(subi).good_resps_cond{condi});
    behavStats.mdWT(subi,condi)= nanmedian(pickupBehav(subi).good_resps_cond{condi});
    behavStats.stdWT(subi,condi)= nanstd(pickupBehav(subi).good_resps_cond{condi});
    behavStats.semWT(subi,condi)= sem(pickupBehav(subi).good_resps_cond{condi});
    behavStats.minWT(subi,condi)= nanmin(pickupBehav(subi).good_resps_cond{condi});
    behavStats.maxWT(subi,condi)= nanmax(pickupBehav(subi).good_resps_cond{condi});
   
    end
end

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

% Grandaverages
GAVGbehav= [];
for condi= 1:5
    
    GAVGbehav.mWT(condi)= nanmean(behavStats.mWT(:,condi),1);
    GAVGbehav.mdWT(condi)= nanmean(behavStats.mdWT(:,condi),1);
    GAVGbehav.stdWT(condi)= nanmean(behavStats.stdWT(:,condi),1);
    GAVGbehav.semWT(condi)= nanmean(behavStats.semWT(:,condi),1);
    
    GAVGbehav.minWT(condi)= nanmean(behavStats.minWT(:,condi),1);
    GAVGbehav.maxWT(condi)= nanmean(behavStats.maxWT(:,condi),1);
    
end

GAVGLogbehav= [];
for condi= 1:5
    
    GAVGLogbehav.mWT(condi)= nanmean(LogBehavStats.mWT(:,condi),1);
    GAVGLogbehav.mdWT(condi)= nanmean(LogBehavStats.mdWT(:,condi),1);
    GAVGLogbehav.stdWT(condi)= nanmean(LogBehavStats.stdWT(:,condi),1);
    GAVGLogbehav.semWT(condi)= nanmean(LogBehavStats.semWT(:,condi),1);
    
    GAVGLogbehav.minWT(condi)= nanmean(LogBehavStats.minWT(:,condi),1);
    GAVGLogbehav.maxWT(condi)= nanmean(LogBehavStats.maxWT(:,condi),1);
    
end

save DescriptiveStats behavStats  LogBehavStats GAVGbehav GAVGLogbehav;

%% Plots
% Histograms for the response times distribution across participants 

% Plots as a function of condition
% linear (?)
a=2; r=2;n=5;
s = a*r.^(0:n-1); %Function rule for Recursive sequence (24 Oct)

% log
z= 1:5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];

hfigure=figure;
for subi= 1: length(RPsubjs)
    
    plot(s,behavStats.mdWT(subi,:),':o','Linewidth',2,...
    'MarkerSize',5);
    hold on;
    
end

% legend(h1,{'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22'},'Location','northwest');

hold on;
% function H=shadedErrorBar_seb(x,y,errBar,lineProps,transparent)
% example: shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'}); 
h2=shadedErrorBar_seb(s,behavStats.mWT(RPsubjs,:),{@mean,@sem},{'-or','LineWidth',3,'MarkerEdgeColor','r',...
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean'},1);

hold on;
h3=shadedErrorBar_seb(s,behavStats.mdWT(RPsubjs,:),{@median,@sem},{'-ob','LineWidth',3,'MarkerEdgeColor','b',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'},1);

        % errorbar(s,GAVGbehav.mWT,GAVGbehav.semWT,'-or','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean');  % RED= mean
        % hold on
        % errorbar(s,GAVGbehav.mdWT,GAVGbehav.semWT,'-sqb','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'); % BLUE= median 

% legend('SEM','mean','SEM','median','Location','northwest');
% legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','SEM','mean','SEM', 'median'},'Location','northwest');
legend({'subj 03','subj 06','subj 07','subj 08','subj 10','subj 13','subj 15','subj 17','subj 18','subj 19','subj 20','subj 21','SEM','mean','SEM', 'median'},'Location','northwest');

set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Average Waiting Times (sec)');
% xlim([0 34]); % added 24 Oct
title(['Mean & Median Waiting Times (N=' num2str(length(RPsubjs)) ')']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for subi= 1: length(RPsubjs)
    
    
    plot(s,LogBehavStats.mdWT(subi,:),':o');
    hold on;
    
end
hold on;
%  legend(h1,{'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22'},'Location','northwest');

hold on;
% function H=shadedErrorBar_seb(x,y,errBar,lineProps,transparent)
% example: shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'}); 
h2=shadedErrorBar_seb(s,LogBehavStats.mWT(RPsubjs,:),{@mean,@sem},{'-or','LineWidth',3,'MarkerEdgeColor','r',...
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean'},1);

hold on;
h3=shadedErrorBar_seb(s,LogBehavStats.mdWT(RPsubjs,:),{@median,@sem},{'-ob','LineWidth',3,'MarkerEdgeColor','b',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'},1);

        % errorbar(s,GAVGbehav.mWT,GAVGbehav.semWT,'-or','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean');  % RED= mean
        % hold on
        % errorbar(s,GAVGbehav.mdWT,GAVGbehav.semWT,'-sqb','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'); % BLUE= median 

% legend('SEM','mean','SEM','median','Location','northwest');
% legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','SEM','mean','SEM', 'median'},'Location','northwest');
legend({'subj 03','subj 06','subj 07','subj 08','subj 10','subj 13','subj 15','subj 17','subj 18','subj 19','subj 20','subj 21','SEM','mean','SEM', 'median'},'Location','northwest');

set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Average Waiting Times (sec)');
% xlim([0 34]); % added 24 Oct
title(['Mean & Median Waiting Times (N=' num2str(length(RPsubjs)) ')']);

%%
figure;
for subi= 1: nSubjs
    
    plot(s,LogBehavStats.mdWT(subi,:),'o');
    hold on;
    
end
hold on;
% plot(s,GbehavStats.mWT,':or','LineWidth',2);
errorbar(s,.mWT,GAVGbehav.semWT,'-or','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean')  % RED= mean
hold on
errorbar(s,GAVGbehav.mdWT,GAVGbehav.semWT,'-sqb','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median') % BLUE= median 

set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Waiting Times (sec)');
legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','mean','median'},'Location','northwest');
xlim([0 34]); % added 24 Oct
title('Mean & Median Waiting Times (N=22)'); 


boxplot(behavStats.mWT);

figure;
boxplot(behavStats.mdWT);
z= 1:5;
set(gca,'xtick',z, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Waiting Times (sec)');
title('Waiting Times per condition across all participants (N=22)'); 

figure;
for condi= 1: 5
    
    plot(behavStats.stdWT(:,condi),'-o');
    hold on;
    
end

% set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
% % set(gca,'ytick',behavStats.stdWT(:,1), 'yticklabel',{'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22'});
% xlabel('Conditions (sec)');
% ylabel('Waiting Times (sec)');
% legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','mean','median'},'Location','northwest');
% xlim([0 34]); % added 24 Oct
% title('Std Waiting Times (N=22)');
