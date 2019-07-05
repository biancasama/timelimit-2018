%% CLUSTER-based permutation tests on ERPs
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on April/May 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: compute the 

    OUTPUT: R{condi}or P{condi} in 

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
    
correlation_folder= [results_Path, '/Correlations']; % it can be also current_subj_folder
    if ~exist(fullfile(correlation_folder)); mkdir(fullfile(correlation_folder)); end;  
    
regression_folder= [results_Path, '/Regressions']; % it can be also current_subj_folder
    if ~exist(fullfile(regression_folder)); mkdir(fullfile(regression_folder)); end;
    
statistics_folder= [results_Path, '/Statistics']; % it can be also current_subj_folder
    if ~exist(fullfile(statistics_folder)); mkdir(fullfile(statistics_folder)); end;

%% 1) Time-series CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data
% load 'R_all_ER1'; load 'R_all_ER2'; 
% load 'R_all_ER1Log'; load 'R_all_ER2Log';
load 'R_all_ER1Log_semiLogX'; load 'R_all_ER2Log_semiLogX'; % correggi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s

RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
nSubjs= length(RPsubjs);
% R_all=R_all';

R_some= {};
for subi=1:nSubjs;
    R_some{subi}= R_all{RPsubjs(subi)}
end

% R_some= R_some';


isequal(R_all{13},R_some{6}); % IT WORKS but you should automatize it

% Run cluster test

                % ROI for RP
%                 chans = {'EEG20','EEG21','EEG29','EEG30','EEG31','EEG39','EEG40'};
                % single-channel
                % chans= {'EEG030'}; 
                chans= {'EEG*'};
                latency= [-2  -0.2]; %-200ms = point of no return;
                nSubjs = 22;
                % Prepare the 'design'. Since Fieldrip does not know the the simple t-test
                % against 0, we have to use the paired t-test (against 0).
                design = zeros(2, 2*nSubjs);
                design(1, 1:nSubjs)= 1;               % condition 1 (data)
                design(1, (nSubjs+1):(2*nSubjs))= 2;        % condition 2 (0, value against which we compare the data).
                design(2, 1:nSubjs) = 1:nSubjs;          % Indicate the subject number (for pairing, alghough not useful here)
                design(2, (nSubjs+1):(2*nSubjs)) = 1:nSubjs;

                % Make *fake data*
                R_contr= {};
                for subi=1:nSubjs;

                    templateER.avg = R_all{subi}.avg; % some or all
                    R_contr{subi}= templateER;
                    R_contr{subi}.avg= zeros(nChans,nTimes);
                end

                cfg = [];
                cfg.method = 'montecarlo';
                cfg.channel = chans;
                cfg.latency= latency;
                cfg.avgoverchan = 'yes';
                % cfg.avgovertime= 'yes';
                cfg.correctm = 'cluster';
                cfg.clusteralpha = 0.05; % lower to 0.01
                cfg.correcttail = 'alpha'; %only for two sided t-tests
                cfg.statistic = 'ft_statfun_depsamplesT';
                % cfg.statistic = 'ft_statfun_indepsamplesF'; % if you want
                % to run the ANOVA
                cfg.design = design;
                cfg.ivar = 1;
                cfg.uvar = 2;
                cfg.numrandomization = 500; % to try put 500 and final put 1000

                % aaron's way
                cfg_neighb=[];
                % cfg_neighb.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
                cfg_neighb.layout= 'eeg_64_NM20884N.lay';
                cfg_neighb.method = 'distance';
                neighbs = ft_prepare_neighbours(cfg_neighb);
                cfg.neighbours = neighbs;

                % STATS  
                stat = ft_timelockstatistics(cfg, R_all{:}, R_contr{:}); %R_some or R_all

                cd(statistics_folder);
%                 save stat_corr_ER1b stat; save stat_corr_ER2b stat; save stat_corr_ER2Logb stat; save stat_corr_ER2Logc stat; save stat_corr_ER2Loge stat;
%                 save stat_corr_ER2Logf stat; save stat_corr_ER2Logg stat; save stat_corr_ER2Logh stat;
                save('stat_corr_ER2Log_wd_22','-v7.3'); % wd= wide: All channels, latency= -3s 0s
                save('stat_corr_ER2Log_nr_22','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_corr_ER2Log_hp_22','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s
                
                save('stat_corr_ER1Log_wd_12','-v7.3'); % wd= wide: All channels, latency= -3s 0s
                save('stat_corr_ER1Log_nr_12','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_corr_ER1Log_hp_12','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s

                stat.posclusters; stat.negclusters;
                
          
                %% Plot significant clusters
                
                stat.stat2= stat.stat.*(stat.posclusterslabelmat == 1);
                cfg_fig= [];
                cfg_fig.parameter= 'stat2';
                cfg_fig.zlim= 'maxabs';
                figure; ft_singleplotER(cfg_fig,stat);
                % 
                % stat.stat_neg= stat.stat.*(stat.negclusterslabelmat == 1);
                % cfg= [];
                % cfg.parameter= 'stat3';
                % cfg.zlim= 'maxabs';
                % figure; ft_singleplotTFR(cfg,stat);
                % % figure; ft_topoplotTFR(cfg,stat);
                
                
                
                cfg = [];
%                 cfg.alpha  = 0.05;
                % cfg.parameter = 'raweffect';
                % cfg.zlim   = [-1e-27 1e-27];
                cfg.layout = 'eeg_64_NM20884N.lay';
                
%                 cfg.saveaspng= 'clusterplot4';  
                ft_clusterplot(cfg,stat);

                figure;
                % define parameters for plotting
                timestep = 0.05;      %(in seconds)
                sampling_rate = 500;
                sample_count = length(stat.time);
                j = [-1:timestep:0];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
                m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
                % get relevant (significant) values
                pos_cluster_pvals = [stat.posclusters(:).prob];

                % In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
                if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.

                pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
                pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

                % First ensure the channels to have the same order in the average and in the statistical output.
                % This might not be the case, because ft_math might shuffle the order
                [i1,i2] = match_str(R_all{1}.label, stat.label);

                % plot
                for k = 1:20;
                    subplot(4,5,k);
                    cfg = [];
                    cfg.xlim=[j(k) j(k+1)];
%                     cfg.zlim = [-5e-14 5e-14];
                    pos_int = zeros(numel(R_all{1}.label),1);
                    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
                    cfg.highlight = 'on';
                    cfg.highlightchannel = find(pos_int);
                    cfg.comment = 'xlim';
                    cfg.commentpos = 'title';
                    cfg.layout = 'eeg_64_NM20884N.lay';
                    ft_topoplotER(cfg, GAVG_R);
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signrank test (move somewhere else?)

                RPSubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
                nRP= length(RPSubjs);
                R_CUBE= zeros(60,2001,nRP);

                for subi=1:RPSubjs;

                    R_CUBE(:,:,subi) = myft_ftavg2dataCube(R_all{subi});

                end

                TAIL= 'both';
                for i=1:nChans

                            for j= 1:nTimes

                [P(i,j)] = signrank(squeeze(R_CUBE(i,j,:)));

                            end

                end

                SR=500;
                t0= 3;
                tAx = [0:2000]./SR - 3;
                figure; plot(tAx, P(28,:), 'o','Linewidth',2);
                title('Correlations between RP average amplitudes from -1s to -0.2s and participants Waiting Times (N=22)');
                Xlabel('Time(s)'); Ylabel('P values');%legend('Correlations p-values','P-value= 0.05','Location','Best');
                hold on; plot([tAx(1) tAx(end)],[0.05 0.05],'r','Linewidth',2);
                txt = 'P-value= 0.05';
                text(-2, 0.07,txt,'FontSize',14);
                plot([0 0],[0 1],'k--','Linewidth',1.5);
                
                %%%%%%%%%%%%%%%%%
                
%                 RPSubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
%                 nRP= length(RPSubjs);
%                 B_CUBE= zeros(60,2001,nRP);
% 
%                 for subi=1:RPSubjs;
% 
%                     B_CUBE(:,:,subi) = myft_ftavg2dataCube(B_all{subi});
% 
%                 end
% 
%                 TAIL= 'both';
%                 for i=1:nChans
% 
%                             for j= 1:nTimes
% 
%                 [P(i,j)] = signrank(squeeze(B_CUBE(i,j,:)));
% 
%                             end
% 
%                 end
% 
%                 SR=500;
%                 t0= 3;
%                 tAx = [0:2000]./SR - 3;
%                 figure; plot(tAx, P(28,:), 'o','Linewidth',2);
%                 title('Correlations between RP average amplitudes from -1s to -0.2s and participants Waiting Times (N=22)');
%                 Xlabel('Time(s)'); Ylabel('P values');%legend('Correlations p-values','P-value= 0.05','Location','Best');
%                 hold on; plot([tAx(1) tAx(end)],[0.05 0.05],'r','Linewidth',2);
%                 txt = 'P-value= 0.05';
%                 text(-2, 0.07,txt,'FontSize',14);
%                 plot([0 0],[0 1],'k--','Linewidth',1.5);

%% 1b) Time-series: REGRESSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data
load 'B_all_ER1'; load 'B_all_ER2'; 
load 'B_all_ER1Log'; load 'B_all_ER2Log';
load 'B_all_ER1Log_semiLogX'; load 'B_all_ER2Log_semiLogX';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s

RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
nSubjs= length(RPsubjs);
B_all=B_all';

B_some= {};
for subi=1:nSubjs;
    B_some{subi}= B_all{RPsubjs(subi)}
end

B_some= B_some';


isequal(B_all{13},B_some{6}); % IT WORKS but you should automatize it

% Run cluster test

                % ROI for RP
%                 chans = {'EEG20','EEG21','EEG29','EEG30','EEG31','EEG39','EEG40'};
                % single-channel
                % chans= {'EEG030'}; 
                chans= {'EEG*'};
                latency= [-2  -0.2]; %-200ms = point of no return;
                nSubjs = 22;
                % Prepare the 'design'. Since Fieldrip does not know the the simple t-test
                % against 0, we have to use the paired t-test (against 0).
                design = zeros(2, 2*nSubjs);
                design(1, 1:nSubjs)= 1;               % condition 1 (data)
                design(1, (nSubjs+1):(2*nSubjs))= 2;        % condition 2 (0, value against which we compare the data).
                design(2, 1:nSubjs) = 1:nSubjs;          % Indicate the subject number (for pairing, alghough not useful here)
                design(2, (nSubjs+1):(2*nSubjs)) = 1:nSubjs;

                % Make *fake data*
                B_contr= {};
                for subi=1:nSubjs;

                    templateER.avg = B_all{subi}.avg;
                    B_contr{subi}= templateER;
                    B_contr{subi}.avg= zeros(nChans,nTimes)
                end

                cfg = [];
                cfg.method = 'montecarlo';
                cfg.channel = chans;
                cfg.latency= latency;
%                 cfg.avgoverchan = 'yes';
                % cfg.avgovertime= 'yes';
                cfg.correctm = 'cluster';
                cfg.clusteralpha = 0.01;
                cfg.correcttail = 'alpha'; %only for two sided t-tests
                cfg.statistic = 'ft_statfun_depsamplesT';
                % cfg.statistic = 'ft_statfun_indepsamplesF'; % if you want
                % to run the ANOVA
                cfg.design = design;
                cfg.ivar = 1;
                cfg.uvar = 2;
                cfg.numrandomization = 500;

                % aaron's way
%                 cfg_neighb=[];
%                 % cfg_neighb.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
%                 cfg_neighb.layout= 'eeg_64_NM20884N.lay';
%                 cfg_neighb.method = 'distance';
%                 neighbs = ft_prepare_neighbours(cfg_neighb);
%                 cfg.neighbours = neighbs;
                
                % elie's way
                cfg_neighb=[];
                % cfg_neighb.layout = '/Users/bt_neurospin/Repos/matlab_internal/turbo_mne/eeg_64_NM20884N.lay';
                cfg_neighb.layout= 'eeg_64_NM20884N.lay';
                cfg_neighb.method = 'distance';
                cfg_neighb.neighbourdist = .18;%.18 %.13
                cfg_neighb.sens= B_all{2}.elec;
                neighbs = ft_prepare_neighbours(cfg_neighb);
                cfg.neighbours = neighbs;

                % STATS  
                stat = ft_timelockstatistics(cfg, B_all{:}, B_contr{:}); %R_all

                cd(statistics_folder); 
%                 save stat_regr_ER1b stat;
%                 save stat_regr_ER2b stat;
%                 save stat_regr_ER2Logb stat;
%                 save stat_regr_ER2Logc stat;
%                 save stat_regr_ER2Loge stat;
%                 save stat_regr_ER2Logf stat;
%                 save stat_regr_ER2Logg stat;
%                 save stat_regr_ER2Logh stat;

                save('stat_regr_ER2Log_wd_22','-v7.3'); % wd= wide: All channels, latency= -3s 0s
                save('stat_regr_ER1Log_nr_22','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_regr_ER1Log_hp_22','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s
                
                save('stat_regr_ER2Log_wd_12','-v7.3'); % wd= wide: All channels, latency= -3s 0s
                save('stat_regr_ER1Log_nr_12','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_regr_ER1Log_hp_12','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s


                stat.posclusters; stat.negclusters;
                
          
                % Plot significant clusters
                
                stat.stat2= stat.stat.*(stat.posclusterslabelmat == 1);
                cfg_fig= [];
                cfg_fig.parameter= 'stat2';
                cfg_fig.zlim= 'maxabs';
                figure; ft_singleplotER(cfg_fig,stat);
                % 
                % stat.stat_neg= stat.stat.*(stat.negclusterslabelmat == 1);
                % cfg= [];
                % cfg.parameter= 'stat3';
                % cfg.zlim= 'maxabs';
                % figure; ft_singleplotTFR(cfg,stat);
                % % figure; ft_topoplotTFR(cfg,stat);

                cfg = [];
                cfg.alpha  = 0.05;
                % cfg.parameter = 'raweffect';
                % cfg.zlim   = [-1e-27 1e-27];
                cfg.layout = 'eeg_64_NM20884N.lay';
%                 cfg.saveaspng= 'clusterplot4';  
                ft_clusterplot(cfg,stat);

                figure;
                % define parameters for plotting
                timestep = 0.05;      %(in seconds)
                sampling_rate = 500;
                sample_count = length(stat.time);
                j = [-1:timestep:0];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
                m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
                % get relevant (significant) values
                pos_cluster_pvals = [stat.posclusters(:).prob];

                % In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
                if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.

                pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
                pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

                % First ensure the channels to have the same order in the average and in the statistical output.
                % This might not be the case, because ft_math might shuffle the order
                [i1,i2] = match_str(R_all{1}.label, stat.label);

                % plot
                for k = 1:20;
                    subplot(4,5,k);
                    cfg = [];
                    cfg.xlim=[j(k) j(k+1)];
%                     cfg.zlim = [-5e-14 5e-14];
                    pos_int = zeros(numel(R_all{1}.label),1);
                    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
                    cfg.highlight = 'on';
                    cfg.highlightchannel = find(pos_int);
                    cfg.comment = 'xlim';
                    cfg.commentpos = 'title';
                    cfg.layout = 'eeg_64_NM20884N.lay';
                    ft_topoplotER(cfg, GAVG_B);
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signrank test (move somewhere else?)

                RPSubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
                nRP= length(RPSubjs);
                R_CUBE= zeros(60,2001,nRP);

                for subi=1:RPSubjs;

                    R_CUBE(:,:,subi) = myft_ftavg2dataCube(R_all{subi});

                end

                TAIL= 'both';
                for i=1:nChans

                            for j= 1:nTimes

                [P(i,j)] = signrank(squeeze(R_CUBE(i,j,:)));

                            end

                end

                SR=500;
                t0= 3;
                tAx = [0:2000]./SR - 3;
                figure; plot(tAx, P(28,:), 'o','Linewidth',2);
                title('Correlations between RP average amplitudes from -1s to -0.2s and participants Waiting Times (N=22)');
                Xlabel('Time(s)'); Ylabel('P values');%legend('Correlations p-values','P-value= 0.05','Location','Best');
                hold on; plot([tAx(1) tAx(end)],[0.05 0.05],'r','Linewidth',2);
                txt = 'P-value= 0.05';
                text(-2, 0.07,txt,'FontSize',14);
                plot([0 0],[0 1],'k--','Linewidth',1.5);
                
                %%%%%%%%%%%%%%%%%
                
%                 RPSubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
%                 nRP= length(RPSubjs);
%                 B_CUBE= zeros(60,2001,nRP);
% 
%                 for subi=1:RPSubjs;
% 
%                     B_CUBE(:,:,subi) = myft_ftavg2dataCube(B_all{subi});
% 
%                 end
% 
%                 TAIL= 'both';
%                 for i=1:nChans
% 
%                             for j= 1:nTimes
% 
%                 [P(i,j)] = signrank(squeeze(B_CUBE(i,j,:)));
% 
%                             end
% 
%                 end
% 
%                 SR=500;
%                 t0= 3;
%                 tAx = [0:2000]./SR - 3;
%                 figure; plot(tAx, P(28,:), 'o','Linewidth',2);
%                 title('Correlations between RP average amplitudes from -1s to -0.2s and participants Waiting Times (N=22)');
%                 Xlabel('Time(s)'); Ylabel('P values');%legend('Correlations p-values','P-value= 0.05','Location','Best');
%                 hold on; plot([tAx(1) tAx(end)],[0.05 0.05],'r','Linewidth',2);
%                 txt = 'P-value= 0.05';
%                 text(-2, 0.07,txt,'FontSize',14);
%                 plot([0 0],[0 1],'k--','Linewidth',1.5);  