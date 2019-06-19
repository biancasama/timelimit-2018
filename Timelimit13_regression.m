%% Regression on time-frequency data
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on April 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the 

    OUTPUT: 

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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load & apply regress/corrcoeff

nChans= 60;
nFreqs= 40;
nTimes= 61;
R= zeros(nChans,nFreqs,nTimes); P= zeros(nChans,nFreqs,nTimes);

powerspectra_folder= [results_Path, '/Powerspect']; % it can be also current_subj_folder
if ~exist(fullfile(powerspectra_folder)); mkdir(fullfile(powerspectra_folder)); end;
cd(powerspectra_folder);

data_Path= powerspectra_folder; % need to fix it in the start up script 

%% Cxorrcoeff 
for subi=1:nSubjs;
    cd(powerspectra_folder);
    load(sprintf('subj%02d_TFR_cond_trl',subi)); %'subj%02d_TFR_bytrial'
    load(sprintf('subj%02d_usefulinfo',subi));
    
    newcond(newcond==Inf) = 32;
    newcond(newcond==2) = NaN;
   
    X= [newcond' ones(length(newcond'),1)];
    Y= TFR_trl.powspctrm;
    if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;
    
%     profile on
    
    R(:)= 0; P(:)=0;
    
    for i=1:nChans
        for k= 1:nFreqs
            for j= 1:nTimes
                
               % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
               [tmpR,tmpP] = corrcoef(Y(:,i,k,j),X(:,1));
               R(i,k,j) = tmpR(1,2);
                
            end
        end
        
        disp(['We are at channel' num2str(i)]);
        
    end
    
    regression_folder= [results_Path, '/Regressions']; % it can be also current_subj_folder
    if ~exist(fullfile(regression_folder)); mkdir(fullfile(regression_folder)); end;
    cd(regression_folder);

    filename= [sprintf('subj%02d_Corr_no2', subi)]; % add one if all trials mixed by condition
    save(filename,'R','P','-v7.3');
    
%     profile off;
%     profile viewer;
    disp(['Subject ' num2str(subi) ' done']);
    
end

%% Grandavg corrcoeff

% make a matrix with all the subjects
load subj22_TFR_oneCorrect
template= TFR_one; %it is subj 22 btw (remember diff in times)
% 
% % Betas_all= {};
% for subi=1:nSubjs;
%     
%     fname_Betas= sprintf('subj%02d_Beta2',subi);
%     pickupBetas(subi) = load(fname_Betas);
%     
%     template.powspctrm = squeeze(pickupBetas(subi).B(1,:,:,:));
%     Betas_all{subi}= template;
%     
% end

R_all= {};
for subi=1:nSubjs;
    
    fname_Corr= sprintf('subj%02d_Corr',subi);
    pickupCorrs(subi) = load(fname_Corr);
    
    template.powspctrm = squeeze(pickupCorrs(subi).R);
    R_all{subi}= template;
    
end

% GAVG_Betas= ft_freqgrandaverage([], Betas_all{:});
GAVG_R= ft_freqgrandaverage([], R_all{:});

cfg=[];
% cfg.baseline= [-2 -1.5];
% cfg.baselinetype= 'relative';
% cfg.zlim= [-3e-25 3e-25];
cfg.xlim= [-2.75 0];
cfg.ylim= [0 30];
% cfg.ylim= [15 30];
cfg.layout = 'eeg_64_NM20884N.lay';
% cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 2;
cfg.showlabels= 'yes';
figure
ft_multiplotTFR(cfg,GAVG_R);

%% STATS 
% chans = {'EEG038','EEG039','EEG040','EEG041','EEG045','EEG046','EEG047','EEG048','EEG049','EEG050','EEG051','EEG053','EEG054','EEG055','EEG056','EEG057','EEG058','EEG059','EEG060'}; %{'EEG028'};
chans= {'EEG*'};
freqs = [14 30];
latency= [-2 -0];
nsubs = 22;
% design = [ones(1,nSubjs) ones(1,nSubjs)*2 ;...
%     1:nSubjs 1:nSubjs];

% Prepare the 'design'. Since Fieldrip does not know the the simple t-test
% against 0, we have to use the paired t-test (against 0).
design = zeros(2, 2*nSubjs);
design(1, 1:nSubjs)= 1;               % condition 1 (data)
design(1, (nSubjs+1):(2*nSubjs))= 2;        % condition 2 (0, value against which we compare the data).
design(2, 1:nSubjs) = 1:nSubjs;          % Indicate the subject number (for pairing, alghough not useful here)
design(2, (nSubjs+1):(2*nSubjs)) = 1:nSubjs;

% Make fake data 

R_contr= {};
for subi=1:nSubjs;
    
    template.powspctrm = R_all{subi}.powspctrm;
    R_contr{subi}= template;
    R_contr{subi}.powspctrm= zeros(nChans,nFreqs,nTimes)
end

cfg = [];
cfg.method = 'montecarlo';
cfg.channel = chans;
cfg.frequency= freqs;
cfg.latency= latency;
% cfg.latency= [-1.1];
% cfg.avgoverchan = 'yes';
% cfg.avgovertime= 'yes';
cfg.frequency = freqs;
cfg.avgoverfreq = 'yes';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.correcttail = 'alpha'; %only for two sided t-tests
cfg.statistic = 'ft_statfun_depsamplesT';
% cfg.statistic = 'ft_statfun_indepsamplesF';
cfg.design = design;
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 10000;

% aaron's way
cfg_neighb=[];
% cfg_neighb.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
cfg_neighb.layout= 'eeg_64_NM20884N.lay';
cfg_neighb.method = 'distance';
neighbs = ft_prepare_neighbours(cfg_neighb);
cfg.neighbours = neighbs;

% % 
stat = ft_freqstatistics(cfg, R_all{:}, R_contr{:});

save stat_try10 stat

stat.posclusters; stat.negclusters;
stat.freq(stat.posclusterslabelmat==1);
stat.time(stat.negclusterslabelmat==1);


stat.stat2= stat.stat.*(stat.posclusterslabelmat == 1);
cfg_fig= [];
cfg_fig.parameter= 'stat2';
cfg_fig.zlim= 'maxabs';
figure; ft_singleplotTFR(cfg_fig,stat);

stat.stat3= stat.stat.*(stat.negclusterslabelmat == 1);
cfg= [];
cfg.parameter= 'stat3';
cfg.zlim= 'maxabs';
figure; ft_singleplotTFR(cfg,stat);
% figure; ft_topoplotTFR(cfg,stat);

cfg = [];
cfg.alpha  = 0.05;
% cfg.parameter = 'raweffect';
% cfg.zlim   = [-1e-27 1e-27];
cfg.layout = 'eeg_64_NM20884N.lay';
ft_clusterplot(cfg,stat);

