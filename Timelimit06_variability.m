%% Compute average & variability within subject
%=========================================================================%
% AUTHORs: Aaron & Yvonne; Bianca modified
% DATE: February/March 2019

%{

    OUTPUT that you need: SF_timecourses_bl (coming from trl + baseline),
    msf and d.


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

%% average 
% make use of DATA_REJ_INTERP
% good_trls = setdiff([1:length(DATA_REJ_INTERP.trial)],allbad_trls); % length of the number of trials kept 
good_trls = setdiff([1:length(DATA_REJ_INTERP.trial)],idx_allbadtrls);
isequal(idx_goodtrls',good_trls)

cond= [TRIALS.cond]; %we put all the conditions in a row
cond(cond==32) = Inf;
un_conds = unique(cond);
newcond= cond(good_trls);

%%
cfg=[];
cfg.trials = good_trls;
% cfg.lpfilter= 'yes';
%%%%%%%%%%%%%%%%%%
%cfg.lpfreq = 30;
%%%%%%%%%%%%%%%%%%
% cfg.lpfreq = 2;
%%%%%%%%%%%%%%%%%%
cfg.demean='yes';
% cfg.baselinewindow = [-0.005 0.005]; % 50 ms but Haggard: 10ms [-5 5]
DATA= ft_preprocessing(cfg,DATA_REJ_INTERP);
% DATA_bl = ft_preprocessing(cfg,DATA_REJ_INTERP); %GOODDATA

avg=[];
for condi = 1:length(un_conds)
    cfg=[];
    cfg.trials= find(newcond == un_conds(condi));
    avg{condi} = ft_timelockanalysis(cfg,DATA);
%     avg{condi} = ft_timelockanalysis(cfg,DATA_bl);
%     mean_avg{condi} = avg{condi}.avg;
end

check= isequal(avg{1}.avg,avg{2}.avg,avg{3}.avg,avg{4}.avg,avg{5}.avg)

%% plot avg 

f1=figure
cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.preproc.lpfilter= 'yes';
cfg.preproc.lpfreq= 35;
cfg.linewidth = 2;
ft_multiplotER(cfg,avg{1},avg{2},avg{3},avg{4},avg{5})

% figure
% cfg= [];
% cfg.layout = 'eeg_64_NM20884N.lay';
% cfg.preproc.lpfilter= 'yes';
% cfg.preproc.lpfreq= 35;
% cfg.channel= [20]
% cfg.linewidth= 2.0;
% ft_singleplotER(cfg,avg{1},avg{2},avg{3},avg{4},avg{5});

%% save stuff

%create folder if it doesn't already exist
[status, msg] = mkdir(current_subj_folder,'Timeseries/Amplitudes/No_Baseline');
cd(fullfile(current_subj_folder,'Timeseries/Amplitudes/No_Baseline'));
% [status, msg] = mkdir(current_subj_folder,'Timeseries/Amplitudes/Baseline_Haggard');
% cd(fullfile(current_subj_folder,'Timeseries/Amplitudes/Baseline_Haggard'));

save timeseries0_EEG avg
% save avg_EEG mean_avg

%% across-trial variability 

for condi = 1:length(un_conds)
    across_stdev{condi} = avg{condi}; % because the fieldtrip function doesn't read in the var
    across_stdev{condi}.avg = sqrt(across_stdev{condi}.var); % we call it avg but it's std
end

%% plot

f2=figure
cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 2;
ft_multiplotER(cfg,across_stdev{1},across_stdev{2},across_stdev{3},across_stdev{4},across_stdev{5})

%% within-trial variability 

within_std = [];
for condi = 1:length(un_conds)
    good_trls= find(newcond == un_conds(condi))
    for i=1:length(good_trls) % length of goodtrls
%         within_std{condi}(i,:) = movstd(DATA.trial{i}(20,:),25); %50
        within_std{condi}(i,:) = movstd(DATA_bl.trial{i}(20,:),25); %50
    end
end

%% plot

for condi = 1:length(un_conds)
    f3(condi)= figure(condi)
    plot(mean(within_std{condi}))
    hold on
end

% change the time axis to put seconds instead of samples

%% Save stuff here

% %create folder if it doesn't already exist
% [status, msg] = mkdir('Timeseries');
% cd(fullfile(current_subj_folder,'/Timeseries'));

% for condi = 1:length(un_conds)
%     avg{condi} = across_stdev{condi}; % because the fieldtrip function doesn't read in the var
%     mean_avg{condi}.avg = avg{condi}.avg;
%     across_stdev{condi}= avg{condi}.var;
% end

% save timeseries_EEG avg
% save avg_EEG mean_avg
save std_EEG across_stdev within_std 

%% %% End (for now)
disp(['END of the script for subj ' int2str(subjnum)])
clear all; close all;
