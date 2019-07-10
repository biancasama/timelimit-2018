% Run cluster test

                % ROI for RP
%                 chans = {'EEG20','EEG21','EEG29','EEG30','EEG31','EEG39','EEG40'};
                % single-channel
                % chans= {'EEG030'}; 
                chans= {'EEG*'};
                latency= [-3  0]; %-200ms = point of no return;
                nSubjs = 12;
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
                stat = ft_timelockstatistics(cfg, R_all{:}, R_contr{:}); %R_some or R_all

                cd(statistics_folder);
%                 save stat_corr_ER1b stat; save stat_corr_ER2b stat; save stat_corr_ER2Logb stat; save stat_corr_ER2Logc stat; save stat_corr_ER2Loge stat;
%                 save stat_corr_ER2Logf stat; save stat_corr_ER2Logg stat; save stat_corr_ER2Logh stat;
                save('stat_corr_ER2Log_wd_22','-v7.3'); % wd= wide: All channels, latency= -3s 0s
                save('stat_corr_ER2Log_nr_22','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_corr_ER1Log_hp_22','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s
                
                save('stat_corr_ER1Log_wd_12','-v7.3'); % wd= wide: All channels, latency= -3s 0s
                save('stat_corr_ER1Log_nr_12','-v7.3'); % nr= narrow: All channels, latency= -2s -0.2s
                save('stat_corr_ER1Log_hp_12','-v7.3'); % hp= hypothesis: ROI channels, latency= -2 -0.2s

                stat.posclusters; stat.negclusters;
                
