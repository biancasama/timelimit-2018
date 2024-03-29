%% Compute Spatial Filter within subject
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

% clear workspace (if needed)

if input('clear all?  (1/0) ... '); clearvars; close all; end
    
% set paths (if needed)
BT_setpath

% choose subj & go to the right folder
BT_getsubj

clear answer LevelAnalysis name numlines prompt subj_folders

%% Load preprocessed files

if subjnum== 1 || subjnum== 17 || subjnum== 18 || subjnum== 19 || subjnum== 20 || subjnum== 21 || subjnum== 22
    load(sprintf('TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp',subjnum))
else
    load(sprintf('TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp',subjnum))
end

%% Clean data more and prepare for Spatial Filter
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

% clear cond un_conds 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg=[];
cfg.trials = good_trls;
% cfg.demean = 'yes';
cfg.lpfilter='yes';
cfg.lpfreq= 30; %notch filter 50
    DATA_CLEANED= ft_preprocessing(cfg,DATA_REJ_INTERP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% per condition 
    DATA_CLcond=[];
    for condi = 1:length(un_conds)

        cfg=[];
        cfg.trials= find(newcond == un_conds(condi));
        DATA_CLcond{condi} = ft_preprocessing(cfg,DATA_REJ_INTERP);

    end

%% Do Spatial Filter here

% Define parameters here (sampling Rate, Timelock, indexes condition)
SR=500;
t0= 3;

if isequal(numel(find(newcond==2)),numel(idx_goodxcond{1}))==1; Ind2 = (newcond'==2); end;
    
if isequal(numel(find(newcond==4)),numel(idx_goodxcond{2}))==1; Ind4 = (newcond'==4); end;
    
if isequal(numel(find(newcond==8)),numel(idx_goodxcond{3}))==1; Ind8 = (newcond'==8); end;
    
if isequal(numel(find(newcond==16)),numel(idx_goodxcond{4}))==1; Ind16 = (newcond'==16); end;
    
if isequal(numel(find(newcond==Inf)),numel(idx_goodxcond{5}))==1; IndInf = (newcond'==Inf); end;

    if isequal(Ind2, Ind4, Ind8, Ind16, IndInf)==0; disp('CORRECT'); else disp('You have done something wrong'); end;
    if isequal(numel(Ind2), numel(Ind4), numel(Ind8), numel(Ind16), numel(IndInf))==1; disp('CORRECT'); else disp('You have done something wrong'); end;

% Ind4 = [TRIALS.cond]' == 4;
% Ind8 = [TRIALS.cond]' == 8;
% Ind16 = [TRIALS.cond]' == 16;
% Ind32 = [TRIALS.cond]' == 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HERE you make the data in 3D structure
DATA_CUBE = myft_ftStruct2dataCube(DATA_CLEANED);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Per condition
DATA_3Dcond=[];
for condi = 1:length(un_conds)
    
    DATA_3Dcond{condi} = myft_ftStruct2dataCube(DATA_CLcond{condi});
    
end


%% Save results of each participant in general folder

readydata_folder= [results_Path,'/Prep_Clean'];
    if ~exist([results_Path, '/Prep_Clean']); mkdir([results_Path, '/Prep_Clean']); end
cd(fullfile(readydata_folder));

filename= [sprintf('subj%02d_DATA_READY', subjnum)];
save(filename, 'DATA_CLEANED','DATA_CLcond');
filename= [sprintf('subj%02d_DATA_CUBES', subjnum)];
save(filename, 'DATA_CUBE', 'DATA_3Dcond');

%% Visualization of ROI around Cz before applying SF

cfg=[];
gavg = ft_timelockanalysis([],DATA_CLEANED);
cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
% already filtered at 30Hz
figure
cfg.linewidth = 2;
ft_multiplotER(cfg,gavg);

% All conditions averaged together
    TrialAvg= mean(DATA_CUBE,3);
    elec_clust= [20,21,29,30,31,39,40]; % cluster of 7 channels selected around Cz
    clustavg= mean(TrialAvg(elec_clust,:));

% Per condition [SOMETHING WRONG]
% for condi = 1:length(un_conds)
%     TrialAvg_cond{condi}= mean(DATA_3Dcond{condi},3);
% %     elec_clust= [20,21,29,30,31,39,40]; % cluster of 7 channels selected around Cz
% %     clustavg_cond= mean(TrialAvg_cond(elec_clust,:));
% end

%% Save results of each participant in general folder

timeseries_folder= [results_Path,'/Timeseries'];
    if ~exist([results_Path, '/Timeseries']); mkdir([results_Path, '/Timeseries']); end
cd(fullfile(timeseries_folder));

filename= [sprintf('subj%02d_Timelockavg', subjnum)];
save(filename, 'gavg','TrialAvg');

%% Plot to see if there is a nice RP

tAx = [0:2000]./SR - 3;
h=figure; plot(tAx, clustavg, 'Linewidth',2);
xlabel('Time (s)');
ylabel('mean Amplitude (\muV)');
title(['Subj ' num2str(subjnum) ', avg channels 20,21,29,30,31,39,40']);
legend('All conditions', 'Location','NorthWest');
% legend('Infinity condition', 'Location','NorthWest');

%  save figure for further comparisons
filename= [sprintf('Readiness_Potential_subj%02d', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);

close; clear h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tAx = [0:2000]./SR - 3;
h=figure; plot(tAx, TrialAvg(40,:), 'Linewidth',2);
% h=figure; plot(tAx, TrialAvg(30,:), 'Linewidth',2);
xlabel('Time (s)');
ylabel('mean Amplitude (\muV)');
title(['Subj ' num2str(subjnum) ', avg Cz']);
legend('All conditions', 'Location','NorthWest');
% legend('Infinity condition', 'Location','NorthWest');

%  save figure for further comparisons
filename= [sprintf('RP_Cz_subj%02d', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);

close; clear h;

%% TOPOPLOT

cfg = [];                            
cfg.xlim = [-1 -0.2]; 
cfg.layout = 'eeg_64_NM20884N.lay';            
h=figure; ft_topoplotER(cfg,gavg); colorbar
title(['Subj ' num2str(subjnum) ', Topography RP']);

%  save figure for further comparisons
filename= [sprintf('Readiness_Potential_topography_subj%02d', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);

close; clear h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define window for Spatial Filter based on the peak amplitude (CRITICAL PART)

% Peak_Lat= cursor_info.Position(1);

% Take window of interest

cfg = [];
cfg.latency = [-0.2 0.2];
Winavg = ft_selectdata(cfg,gavg);

elec_clust= [20,21,29,30,31,39,40];
selclustavg= mean(Winavg.avg(elec_clust,:));
% Find minimum 
[Amp,Lat]= min(selclustavg);
Peak_Lat = Winavg.time(Lat);

h=figure; plot(Winavg.time, selclustavg, 'Linewidth',2);
% click with cursor first 

filename= [sprintf('subj%02d_peakLat', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);

close; clear h;

Peak_Win= [Peak_Lat-0.05 Peak_Lat+0]; %[Peak_Lat-0.05 Peak_Lat+0.05];

% use the time difference to get the spatial filter
tWin= {Peak_Win-1, Peak_Win}; % define tWin based on grand average for that subject

%% Quantitavely assess if there is a RP (microVolt decrease method)
% Alternative method to Slope analysis

M1=mean(clustavg(seconds2samples(tWin{1},t0,SR)));
M2=mean(clustavg(seconds2samples(tWin{2},t0,SR)));
if M2-M1 < -1e-6 % 1 microVolt
    display('OK')
else
    disp('NOT OK')
end

RP_test= M2-M1;

%% SPATIAL FILTER here
% We apply on all the conditions to be fair 
nTrials = size(DATA_CUBE,3);
IndAll = true(nTrials,1);
[trl,msf] = ems_ncond(DATA_CUBE,IndAll,[], @spf_basic_2time_diff,tWin,t0,SR);

%% using the peak to get the spatial filter
% Explain here better which method it is

% d = mean(mean(DATA_CUBE(:,seconds2samples(tWin{2},t0,SR),:),3),2);
% d = d ./ norm(d);
% for i=1:162
%     RESULT(i,:) = -(DATA_CUBE(:,:,i)' * d)';
% end
% plot(tAx,mean(RESULT))
% cfg.layout = 'eeg_64_NM20884N.lay';
% figure;topoplot(cfg,d);
% colorbar

%% Apply baseline peak-based
% blRange = [0 0.050]; % RP peak
blRange= Peak_Win;
% blRange = [0.225 0.275]; % somato-motor potential
SF_timecourses_bl = baseline_correct(trl,SR,3,blRange);

%% Visualization of the time course after the Spatial Filter

h=figure; plot(tAx,mean(SF_timecourses_bl), 'Linewidth',2);
xlabel('Time (s)');
ylabel('mean Amplitude (\muV)');
title(['Subj ' num2str(subjnum) ', Spatial filter, baseline [ ' num2str(tWin{2}) ' ]']);
legend('All conditions', 'Location','best');

% save figure for further comparisons
filename= [sprintf('Spatial_filter_subj%02d', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);
close; clear h;

%% Visualization of the topography after the Spatial Filter

meanSF = mean(msf,2); 
cfg.layout = 'eeg_64_NM20884N.lay';
% cfg.xlim = [-1 -0.2];
h=figure;topoplot(cfg,meanSF);
title(['Subj ' num2str(subjnum) ', Topography Spatial filter, baseline [ ' num2str(tWin{2}) ' ]']);
colorbar

% save figure for further comparisons
filename= [sprintf('Spatial_filter_topography_subj%02d', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);
close; clear h;

% comparison with before SF, sSubj 22ame window

cfg = [];                            
cfg.xlim= [tWin{1}(1) tWin{2}(2)];filename= [sprintf('subj%02d_SF_results_newBl', subjnum)];
save(filename, 'Ind*', 'tWin', 'nTrials', 'trl', 'msf', 'blRange', 'SF_timecourses_bl', 'tAx', 'subjnum', 'RP_test');
cfg.layout = 'eeg_64_NM20884N.lay';            
h=figure; ft_topoplotER(cfg,gavg); colorbar
title(['Subj ' num2str(subjnum) ', Topography RP']);

%  save figure for further comparisons
filename= [sprintf('Readiness_Potential_top_comparison_subj%02d', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);

%% Save results for each participant in their own folder 

% if input('Save SPATIAL FILTER results? ')
%     timeseries_folder= [current_subj_folder,'/Timeseries'];
%     if ~exist([timeseries_folder, '/SpatialFilter']); mkdir([timeseries_folder, '/SpatialFilter']); end
%     cd(fullfile(timeseries_folder,'/SpatialFilter'));
% 
%     save SF_results_newBl Ind* tWin nTrials trl msf blRange SF_timecourses_bl tAx subjnum
% end

%% Save results of each participant in general folder

spatialfilter_folder= [results_Path,'/SpatialFilter'];
    if ~exist([results_Path, '/SpatialFilter']); mkdir([results_Path, '/SpatialFilter']); end
cd(fullfile(spatialfilter_folder));

filename= [sprintf('subj%02d_SF_results_newBl', subjnum)];
save(filename, 'Ind*', 'tWin', 'nTrials', 'trl', 'msf', 'blRange', 'SF_timecourses_bl', 'tAx', 'subjnum', 'RP_test');

%% Tell me what I have done so far
disp(['subj ' num2str(subjnum) ' done']);

close all;
clearvars;

%% End of the script
