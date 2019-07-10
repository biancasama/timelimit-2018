%% Super cluster files 
% Load in useful data

load 'B_all_ER2Log';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s

% Run cluster test

                % ROI for RP
%                 chans = {'EEG20','EEG21','EEG29','EEG30','EEG31','EEG39','EEG40'};
                % single-channel
                % chans= {'EEG030'}; 
                chans= {'EEG*'};
                latency= [-3  0]; %-200ms = point of no return;
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

                    templateER.avg = B_some{subi}.avg;
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
                stat = ft_timelockstatistics(cfg, B_some{:}, B_contr{:}); %R_all

                cd(statistics_folder);
%                 save stat_regr_ER1b stat;
%                 save stat_regr_ER2b stat;
%                 save stat_regr_ER2Logb stat;
%                 save stat_regr_ER2Logc stat;
%                 save stat_regr_ER2Loge stat;
%                 save stat_regr_ER2Logf stat;
%                 save stat_regr_ER2Logg stat;
%                 save stat_regr_ER2Logh stat;

                save('stat_regr_ER1Log_wd_22','-v7.3'); % wd= wide: All channels, latency= -3s 0s
                save('stat_regr_ER1Log_nr_22','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_regr_ER1Log_hp_22','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s
                
                


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
                latency= [-3  0]; %-200ms = point of no return;
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

                    templateER.avg = B_some{subi}.avg;
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
                stat = ft_timelockstatistics(cfg, B_some{:}, B_contr{:});
                
                