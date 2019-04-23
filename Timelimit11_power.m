%% Compute power spectrum within subject
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
            cfg.latency = [-1 -.2];
            DATA_CLEAN= ft_selectdata(cfg,DATA_REJ_INTERP);
  
    % DATA_bl = ft_preprocessing(cfg,DATA_REJ_INTERP); % In case you apply
    % baseline 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                freq=[];
                for condi = 1:length(un_conds)
                    
                    cfg=[];
                    cfg.trials= find(newcond == un_conds(condi));
                    cfg.method= 'mtmfft';
                    cfg.taper= 'hanning';
                    cfg.foilim= [1 40];
                    cfg.pad = 2;
                    %             cfg.width= 7;
                    %             cfg.gwidth= 3;
                    freq{condi} = ft_freqanalysis(cfg,DATA_CLEAN);
                    
                    % Save averages power spectra
                    
                    % Create the folder if it doesn't exist already.
                   disp('Save SPECTRAL AVERAGES results? RISK OF OVERWRITING... ')
                        
                        powerspectra_folder= [results_Path, '/Powerspect']; % it can be also current_subj_folder
                        if ~exist(fullfile(powerspectra_folder)); mkdir(fullfile(powerspectra_folder)); end;
                        cd(powerspectra_folder);
                        
                        filename= [sprintf('subj%02d_RSp_freq', subi)];
                        save(filename,'freq','-v7.3');
            
                    
                end
                
end

cd([results_Path, '/Powerspect']);
disp('END of (TIME)FREQ ANALYSIS by condition per subject');

%% Grandaverage
nSubjs= 22;
%~~~ Other
OKsubjs= [3 6 7 8 10 13 15 17 18 19 20 21];

cd(pwd);
for subi=1:nSubjs;
    
    fname_freq= sprintf('subj%02d_RSp_freq',subi);
    pickupFreq(subi) = load(fname_freq);
   
end


%
freqmatrix=[];
for subi= 1:nSubjs; %nGoodSubjects
    
    for k= 1:5
        freqmatrix{subi,k}= pickupFreq(subi).freq{k}; %pickupSub(i).avg_EEG{k}
    end
end

save freqmatrix freqmatrix

Grand_Freq=[]; 
for k= 1:5
    
%     cfg=[];
%   
%     cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
    Grand_Freq{k}= ft_freqgrandaverage([],freqmatrix{:,k});
%     Grand_Freq_ROI{k}= mean(Grand_Freq{k}.powspctrm,1);
    
end

save Grand_Freq Grand_Freq

cfg=[];
cfg.layout = 'eeg_64_NM20884N.lay';
cfg.linewidth = 2;
cfg.showlabels= 'yes';
figure
ft_multiplotER(cfg,Grand_Freq{1},Grand_Freq{2},Grand_Freq{3},Grand_Freq{4},Grand_Freq{5});

%%
powermatrix_premov= struct('alpha', [], 'beta', []); 
mean_premov_freq= struct('alpha', [], 'beta', []); 

for i= 1: nSubjs; %nGoodSubjects or nSubjects
    
    for k= 1:5
        
        cfg = [];
        cfg.channel = {'EEG038','EEG039','EEG040','EEG041','EEG045','EEG046','EEG047','EEG048','EEG049','EEG050','EEG051','EEG053','EEG054','EEG055','EEG056','EEG057','EEG058','EEG059','EEG060'};
        cfg.frequency= [7.5 12.5];
%         powermatrix_premov.alpha{i,k} = ft_selectdata(cfg, freqmatrix{i,k});
%         powermatrix_premov.alpha{i,k}= mean(powermatrix_premov.alpha{i,k}.powspctrm,1);
        mean_premov_freq.alpha(i, k) = mean(powermatrix_premov.alpha{i,k});
%         (i, k) = mean(datamatrix_premov.alpha{i,k}.powspctrm);
%         mean_premov_freq.alpha
%         datamatrix_premov.ch30{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
%         mean_premov_amp.ch30(i, k) = mean(datamatrix_premov.ch30{i,k}.avg);
    end
end

%%
% Alpha= 0.05,  tail 'left'
TAIL= 'right'; %'both','left','right'.
% Condition n.1 (2sec) vs all the others
[p1,h1,stats] =signrank(mean_premov_freq.alpha(:, 1),mean_premov_freq.alpha(:, 5),'tail',TAIL) % vs Inf
[p2,h2,stats] =signrank(mean_premov_freq.alpha(:, 1),mean_premov_freq.alpha(:, 4),'tail',TAIL) % vs 16sec
[p3,h3,stats] =signrank(mean_premov_freq.alpha(:, 1),mean_premov_freq.alpha(:, 3),'tail',TAIL) % vs 8sec
[p4,h4,stats] =signrank(mean_premov_freq.alpha(:, 1),mean_premov_freq.alpha(:, 2),'tail',TAIL) % vs 4sec

% Condition n.2 (4sec) vs all the others
[p5,h5,stats] =signrank(mean_premov_freq.alpha(:, 2),mean_premov_freq.alpha(:, 5),'tail',TAIL) % vs Inf
[p6,h6,stats] =signrank(mean_premov_freq.alpha(:, 2),mean_premov_freq.alpha(:, 4),'tail',TAIL) % vs 16sec
[p7,h7,stats] =signrank(mean_premov_freq.alpha(:, 2),mean_premov_freq.alpha(:, 3),'tail',TAIL) % vs 8sec

% Condition n.3 (8sec) vs all the others
[p8,h8,stats8] =signrank(mean_premov_freq.alpha(:, 3),mean_premov_freq.alpha(:, 5),'tail',TAIL) % vs Inf
[p9,h9,stats] =signrank(mean_premov_freq.alpha(:, 3),mean_premov_freq.alpha(:, 4),'tail',TAIL) % vs 16sec

% Condition n.4 (16sec) vs 5
[p10,h10,stats] =signrank(mean_premov_freq.alpha(:, 4),mean_premov_freq.alpha(:, 5),'tail',TAIL) % vs Inf

% FDR
pvals_SF_05= [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10];
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals_SF_05,0.05,'pdep','yes');

%%
GROUP= {'2 sec','4 sec', '8 sec','16 sec','Inf'};

% Alpha range 
[P,ANOVATAB,STATS] = anova1(mean_premov_freq.alpha,GROUP);
[c,m,h,nms] = multcompare(STATS);

[P,ANOVATAB,STATS] = kruskalwallis(mean_premov_freq.alpha,GROUP);
[c,m,h,nms] = multcompare(STATS);

%% cluster stats

% use freqmatrix
chans = {'EEG038','EEG039','EEG040','EEG041','EEG045','EEG046','EEG047','EEG048','EEG049','EEG050','EEG051','EEG053','EEG054','EEG055','EEG056','EEG057','EEG058','EEG059','EEG060'};
freqs = [7 14];
nsubs = 22;
design = [ones(1,nsubs) ones(1,nsubs)*2; ...
    1:nsubs 1:nsubs];
mat = freqmatrix';

cfg = [];
cfg.method = 'montecarlo';
cfg.channel = chans;
cfg.avgoverchan = 'yes';
cfg.frequency = freqs;
cfg.avgoverfreq = 'yes';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.correcttail = 'alpha';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.design = design;
cfg.ivar = 1;
cfg.uvar = 2;
cfg.numrandomization = 10000;
stat = ft_freqstatistics(cfg, mat{1,:}, mat{4,:})




