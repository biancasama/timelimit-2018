%% STATS 
%% Cluster-based permutation tests
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

%%
statistics_folder= [results_Path, '/Statistics']; % it can be also current_subj_folder
    if ~exist(fullfile(statistics_folder)); mkdir(fullfile(statistics_folder)); end;

%% 1) Time-series: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data
load 'R_all_ER1'; load 'R_all_ER2'; load 'R_all_ER2Log';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];
nSubjs= length(RPsubjs);
R_all=R_all';

R_some= {};
for subi=1:nSubjs;
    R_some{subi}= R_all{RPsubjs(subi)}
end

R_some= R_some';


isequal(R_all{13},R_some{6}); % IT WORKS but you should automatize it

% Run cluster test

                % ROI for RP
                % chans = {'EEG20','EEG21','EEG29','EEG30','EEG31','EEG39','EEG40'};
                % single-channel
                % chans= {'EEG030'}; 
                chans= {'EEG*'};
                latency= [-3  0]; %-200ms = point of no return;
                nSubs = 12;
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

                    templateER.avg = R_some{subi}.avg;
                    R_contr{subi}= templateER;
                    R_contr{subi}.avg= zeros(nChans,nTimes)
                end

                cfg = [];
                cfg.method = 'montecarlo';
                cfg.channel = chans;
                cfg.latency= latency;
%                 cfg.avgoverchan = 'yes';
                % cfg.avgovertime= 'yes';
                cfg.correctm = 'cluster';
                cfg.clusteralpha = 0.05;
                cfg.correcttail = 'alpha'; %only for two sided t-tests
                cfg.statistic = 'ft_statfun_depsamplesT';
                % cfg.statistic = 'ft_statfun_indepsamplesF'; % if you want
                % to run the ANOVA
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

                % STATS  
                stat = ft_timelockstatistics(cfg, R_some{:}, R_contr{:}); %R_all

                cd(statistics_folder);
                save stat_corr_ER1b stat; save stat_corr_ER2b stat; save stat_corr_ER2Logb stat; save stat_corr_ER2Logc stat; save stat_corr_ER2Loge stat;
                save stat_corr_ER2Logf stat; save stat_corr_ER2Logg stat; save stat_corr_ER2Logh stat;

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
                cfg.saveaspng= 'clusterplot4';  
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
                
%% 1) Time-frequency: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load in useful data
load 'R_all_TFR1'; load 'R_all_TFR2'; load 'R_all_TFR3';
uconds=  [2     4     8    16   Inf];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run cluster test

%% STATS 
% chans = {'EEG038','EEG039','EEG040','EEG041','EEG045','EEG046','EEG047','EEG048','EEG049','EEG050','EEG051','EEG053','EEG054','EEG055','EEG056','EEG057','EEG058','EEG059','EEG060'}; %{'EEG028'};
                chans= {'EEG*'};
%                 freqs = [14 30];
%                 latency= [-2 -0];
                nsubs = 22;
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

                %     template.powspctrm = R_all{subi}.powspctrm;
                %     R_contr{subi}= template;
                %     R_contr{subi}.powspctrm= zeros(nChans,nFreqs,nTimes)

                    for condi = 1:length(uconds)

                        template{condi}.powspctrm = R_all{subi,condi}.powspctrm;
                        R_contr{subi,condi}= template{condi};
                        R_contr{subi,condi}.powspctrm= zeros(nChans,nFreqs,nTimes);

                    end

                end


                cfg = [];
                cfg.method = 'montecarlo';
                cfg.channel = chans;
                cfg.frequency= 'all';
                cfg.latency= 'all';
%                 % cfg.latency= [-1.1];
%                 % cfg.avgoverchan = 'yes';
%                 % cfg.avgovertime= 'yes';
%                 cfg.frequency = freqs;
%                 cfg.avgoverfreq = 'yes';
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
                
                stat= {};

                for condi = 1:length(uconds)

                    stat{condi} = ft_freqstatistics(cfg, R_all{:,condi}, R_contr{:,condi});

                end

                save stat_corr_cond stat

                stat.posclusters; stat.negclusters;
                % stat.freq(stat.posclusterslabelmat==1);
                % stat.time(stat.negclusterslabelmat==1);

                neg_cluster_pvals = [stat{2}.negclusters(:).prob];
                neg_signif_clust = find(neg_cluster_pvals < stat{2}.cfg.alpha);
                neg = ismember(stat{2}.negclusterslabelmat, neg_signif_clust);
                stat{2}.time(neg(29,1,:));


                % stat.stat2= stat.stat.*(stat.posclusterslabelmat == 1);
                % cfg_fig= [];
                % cfg_fig.parameter= 'stat2';
                % cfg_fig.zlim= 'maxabs';
                % figure; ft_singleplotTFR(cfg_fig,stat);
                % 
                stat{2}.stat3= stat{2}.stat*(stat.negclusterslabelmat == 1);
                cfg= [];
                cfg.parameter= 'stat2';
                cfg.zlim= 'maxabs';
                figure; ft_singleplotER(cfg,stat{2});
                % figure; ft_topoplotTFR(cfg,stat);

                cfg = [];
                cfg.alpha  = 0.05;
                % cfg.parameter = 'raweffect';
                cfg.zlim   = [-4.5 4.5];
                cfg.layout = 'eeg_64_NM20884N.lay';
                ft_clusterplot(cfg,stat{2});

