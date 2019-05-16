%% Compute time-frequency within subject
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on April 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the 

    OUTPUT: freq{condi} in 

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
    %             cfg.latency = [-1 -.2];
%     DATA_CLEAN= ft_selectdata(cfg,DATA_REJ_INTERP);
    DATA_CLEAN= ft_preprocessing(cfg,DATA_REJ_INTERP);

    % DATA_bl = ft_preprocessing(cfg,DATA_REJ_INTERP); % In case you apply
    % baseline
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TFR_cond={};
    for condi = 1:length(un_conds)
        
        cfg=[];
        cfg.trials= find(newcond == un_conds(condi));
%         cfg.keeptrials  = 'yes';
        cfg.output= 'pow';
        cfg.method= 'mtmconvol'; %wavelet
        cfg.taper= 'dpss';
        cfg.foi= 1:1:40;
        cfg.toi= -3: 0.05:0;
        cfg.pad = 'nextpow2';
        %                     cfg.width= 7;
        %                     cfg.gwidth= 3;
        cfg.t_ftimwin= ones(length(cfg.foi))*.5;
        cfg.tapsmofrq= linspace(2,6,length(cfg.foi));
        TFR_cond{condi} = ft_freqanalysis(cfg,DATA_CLEAN);
        
    end
%     TFR_one = ft_freqanalysis(cfg,DATA_CLEAN);
    
    % Save averages power spectra
    
    % Create the folder if it doesn't exist already.
%     disp('Save TIMEFREQ AVERAGES results? RISK OF OVERWRITING... ');
   
    powerspectra_folder= [results_Path, '/Powerspect']; % it can be also current_subj_folder
    if ~exist(fullfile(powerspectra_folder)); mkdir(fullfile(powerspectra_folder)); end;
    cd(powerspectra_folder);
    
    filename= [sprintf('subj%02d_TFR_cond', subi)]; % add one if all trials mixed by condition
    save(filename,'TFR_cond','-v7.3');
    
    disp(['Subject ' num2str(subi) ' done']);
    
end

% filename= [sprintf('subj%02d_usefulinfo', subi)]; % add one if all trials mixed by condition
% save(filename,'newcond','-v7.3');
% disp(['Subject ' num2str(subi) ' done']);
% 
% end

cd([results_Path, '/Powerspect']);
disp('END of (TIME)FREQ ANALYSIS by condition per subject');

%% Grandaverage
nSubjs= 22;
%~~~ Other
OKsubjs= [3 6 7 8 10 13 15 17 18 19 20 21];

cd(pwd);
for subi=1:nSubjs;
    
    fname_TFR= sprintf('subj%02d_TFR_cond',subi); %TFR_one
    pickupTFR(subi) = load(fname_TFR);
   
end


%
TFRmatrix=[];
for subi= 1:nSubjs; %nGoodSubjects
    
    for k= 1:5
        TFRmatrix{subi,k}= pickupTFR(subi).TFR_cond{k}; %pickupSub(i).avg_EEG{k}
    end
end

save TFRmatrix TFRmatrix

Grand_TFR=[]; 
for k= 1:5
    
%     cfg=[];
%   
%     cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
    Grand_TFR{k}= ft_freqgrandaverage([],TFRmatrix{:,k});
%     Grand_Freq_ROI{k}= mean(Grand_Freq{k}.powspctrm,1);
    
end

save Grand_TFR Grand_TFR

cfg=[];
cfg.baseline= [-2 -1.5];
cfg.baselinetype= 'relative';
% cfg.zlim= [-3e-25 3e-25];
cfg.xlim= [-1 -0];
% cfg.ylim= [8 30];
% cfg.ylim= [15 30];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 2;
cfg.showlabels= 'yes';
figure
ft_multiplotTFR(cfg,Grand_TFR);%Grand_TFR{2},Grand_TFR{3},Grand_TFR{4},Grand_TFR{5}

cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 2;
cfg.showlabels= 'yes';
figure
ft_multiplotER(cfg,Grand_TFR{1},Grand_TFR{2},Grand_TFR{3},Grand_TFR{4},Grand_TFR{5});

% ft_topoplotTFR(cfg,Grand_TFR{2});
%%
datamatrix_premov= struct('alpha', []); 
mean_premov_TFR= struct('alpha', []); 

for i= 1: nSubjs; %nGoodSubjects or nSubjects
    
    for k= 1:5
        
        cfg = [];
        cfg.channel = {'EEG028'};
        cfg.frequency= [8 30];
        datamatrix_premov.alpha{i,k} = ft_selectdata(cfg, TFRmatrix{i,k});
        mean_premov_TFR.alpha(i, k) = mean(mean(datamatrix_premov.alpha{i,k}.powspctrm));
        
    end
end

%%
% Alpha= 0.05,  tail 'left'
TAIL= 'both'; %'both','left','right'.
% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(mean_premov_TFR.alpha(:, 1),mean_premov_TFR.alpha(:, 5),'tail',TAIL) % vs Inf
[p2,h2,stats] =signrank(mean_premov_TFR.alpha(:, 1),mean_premov_TFR.alpha(:, 4),'tail',TAIL) % vs 16sec
[p3,h3,stats] =signrank(mean_premov_TFR.alpha(:, 1),mean_premov_TFR.alpha(:, 3),'tail',TAIL) % vs 8sec
[p4,h4,stats] =signrank(mean_premov_TFR.alpha(:, 1),mean_premov_TFR.alpha(:, 2),'tail',TAIL) % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(mean_premov_TFR.alpha(:, 2),mean_premov_TFR.alpha(:, 5),'tail',TAIL) % vs Inf
[p6,h6,stats] =signrank(mean_premov_TFR.alpha(:, 2),mean_premov_TFR.alpha(:, 4),'tail',TAIL) % vs 16sec
[p7,h7,stats] =signrank(mean_premov_TFR.alpha(:, 2),mean_premov_TFR.alpha(:, 3),'tail',TAIL) % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(mean_premov_TFR.alpha(:, 3),mean_premov_TFR.alpha(:, 5),'tail',TAIL) % vs Inf
[p9,h9,stats] =signrank(mean_premov_TFR.alpha(:, 3),mean_premov_TFR.alpha(:, 4),'tail',TAIL) % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(mean_premov_TFR.alpha(:, 4),mean_premov_TFR.alpha(:, 5),'tail',TAIL) % vs Inf

% FDR
pvals_SF_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adBetas_allj_p]=fdr_bh(pvals_SF_05,0.05,'pdep','yes');

%%
GROUP= {'2 sec','4 sec', '8 sec','16 sec','Inf'};

% Alpha range 
[P,ANOVATAB,STATS] = anova1(mean_premov_TFR.alpha,GROUP);
[c,m,h,nms] = multcompare(STATS);

[P,ANOVATAB,STATS] = kruskalwallis(mean_premov_TFR.alpha,GROUP);
[c,m,h,nms] = multcompare(STATS);

%% cluster stats

% use freqmatrix
% chans = {'EEG038','EEG039','EEG040','EEG041','EEG045','EEG046','EEG047','EEG048','EEG049','EEG050','EEG051','EEG053','EEG054','EEG055','EEG056','EEG057','EEG058','EEG059','EEG060'}; %{'EEG028'};
chans= {'EEG*'};
freqs = [1 40];
latency= [-3 -1];
nsubs = 22;
% design = [ones(1,nsubs) ones(1,nsubs)*2 ; ...
%     1:nsubs 1:nsubs];
design = [ones(1,nsubs) ones(1,nsubs)*2 ones(1,nsubs)*3 ones(1,nsubs)*4 ones(1,nsubs)*5; ...
    1:nsubs 1:nsubs 1:nsubs 1:nsubs 1:nsubs];
% TFRall = TFRmatrix';

cfg = [];
cfg.method = 'montecarlo';
cfg.channel = chans;
cfg.latency= latency;
cfg.latency= 'all';
cfg.avgoverchan = 'yes';
% cfg.avgovertime= 'yes';
% cfg.frequency = freqs;
% cfg.avgoverfreq = 'yes';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.correcttail = 'alpha'; %only for two sided t-tests
cfg.statistic = 'ft_statfun_depsamplesT';
% cfg.statistic = 'ft_statfun_indepsamplesF';
cfg.design = design;
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 10000;

% avgall= avgmatrix';
% 
% cfg_neighb= [];
% cfg_neighb.method= 'distance';
% cfg_neighb.neighbourdist= 2;
% cfg.channel= {'EEG*'};
% % cfg_neighb.template= 'biosemi64_neighb.mat';
% cfg_neighb.layout = 'eeg_64_NM20884N.lay';
% cfg_neighb.elec= TFRall{1}.elec;
% neighbours= ft_prepare_neighbours(cfg_neighb, avgall{1});
% cfg.neighbours = neighbours;

% cfg_neighb= [];
%     %cfg.layout = '/home/aschurge/matlab/internal/myft/eeg_64_generic.lay';%eeg_64_NM20884N.lay
%     cfg.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
%     cfg.method = 'distance';
%     neighbs = ft_prepare_neighbours(cfg);
%     cfg.neighbours = neighbs;
%     cfg.method = 'nearest';

% aaron's way
cfg_neighb=[];
%cfg.layout = '/home/aschurge/matlab/internal/myft/eeg_64_generic.lay';%eeg_64_NM20884N.lay
cfg_neighb.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
cfg_neighb.method = 'distance';
neighbs = ft_prepare_neighbours(cfg_neighb);
cfg.neighbours = neighbs;

% % 
stat = ft_freqstatistics(cfg, TFRall{1,:},TFRall{2,:},TFRall{3,:},TFRall{4,:},TFRall{5,:});
% 
% stat = ft_freqstatistics(cfg, TFRall{1,:}, TFRall{5,:});
% stat = ft_freqstatistics(cfg, Grand_TFR{1}, Grand_TFR{5});
% 
stat.posclusters; stat.negclusters;
stat.time(stat.posclusterslabelmat==1);
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

%% 


freq2sec= Grand_TFR{1};
freqInf= Grand_TFR{5};

freq2sec= TFRall{1,:};
freqInf= TFRall{5,:};

% 
% cfg.latency= latency;
% cfg.frequency = freqs;
freq2sec= ft_selectdata([],freq2sec);
freqInf= ft_selectdata([],freqInf);

freq2sec= ft_freqdescriptives([],freq2sec);
freqInf= ft_freqdescriptives([],freqInf);

stat.raweffect= freq2sec.powspctrm - freqInf.powspctrm;

% cfg= [];
cfg.alpha = 0.025;
cfg.channel = 'all';
cfg.parameter = 'raweffect';
% cfg.xlim= [-2 -0.2];
% cfg.ylim= [8 30];
% % cfg.zlim = 
cfg.layout = 'eeg_64_NM20884N.lay';
ft_clusterplot([],stat);
