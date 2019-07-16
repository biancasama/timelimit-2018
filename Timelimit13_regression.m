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

%% 1a) Behavioral: REGRESSION between response times (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Response times (alias Waiting Times): containing alredy log
    % transformed data
    
    load(pickupBehav);

    % X= Conditions (alias Time limits): do Log transformation in following
    % loop
    for subi=1:nSubjs;

        cd(powerspectra_folder);
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute regression here:

cd(regression_folder);
B= zeros(nSubjs,1); R= zeros(nSubjs,1); stats= zeros(nSubjs,1);

    for subi=1:nSubjs;

        % Just re-converting Inf values into an integer number 
        pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32; 
        pickupCond(subi).Logcond= log(pickupCond(subi).newcond);
    
        X= [pickupCond(subi).Logcond' ones(length(pickupCond(subi).Logcond'),1)];
        Y= pickupBehav(subi).LogRESPS';
        if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

    %     profile on
% 
%         B(:)= 0; R(:)=0; stats(:)=0;

                   [B,~,R,~,stats]= regress(Y,X);
%                    [tmpR,tmpP] = corrcoef(Y(:),X(:,1));
%                    B= tmpB(1,2); R{:}= tmpR(1,:);

                   filename= [sprintf('subj%02d_Regr_BhvLog', subi)]; % add one if all trials mixed by condition
                   save(filename,'B','R', 'stats','-v7.3');

    %                profile off;
    %                profile viewer;
                   disp(['Subject ' num2str(subi) ' done']);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reload obtained regressions values in a matrix:

    for subi=1:nSubjs;

        fname_Regr= sprintf('subj%02d_Regr_BhvLog',subi);
        pickupRegr(subi) = load(fname_Regr);

    end
    
    save(filename,'pickupRegr_BehavLog', 'pickupRegr','-v7.3');

% Plot regressions:

% ~~~~~
%% 2a) Time-series: REGRESSION between RP amplitudes (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s
% B= zeros(2,nChans,nTimes); R= zeros(nChans,nTimes); stats= zeros(nChans,nTimes);

% clear data_Path; 
% % set new data_Path
% data_Path= timeseries_folder; % need to fix it in the start up script
% 

% Load in useful data

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    load('pickupER');
    
%     for subi=1:nSubjs;
% 
%         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
%         pickupER(subi) = load(fname_ER);
% 
%     end
%     save(filename,'pickupER','pickupER','-v7.3'); % For variables larger than 2GB use MAT-file version 7.3 or later. 
% 

    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    load('pickupCond');
    
%     for subi=1:nSubjs;
%         
%         fname_Cond= sprintf('subj%02d_usefulinfo',subi);
%         pickupCond(subi) = load(fname_Cond);
%         
%     end
%      save 'pickupCond' pickupCond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute regressions here:

cd(regression_folder);

        for subi=1:nSubjs;

            % Just re-converting Inf values into an integer number
            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;
             
            X= [pickupCond(subi).RESPTIMES' ones(length(pickupCond(subi).RESPTIMES'),1)];
            Y= pickupER(subi).avg_trl.trial; %.avg or .trial?
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

           clear B R stats

            for i=1:nChans

                    for j= 1:nTimes

                       [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X);

                    end


                disp(['We are at channel' num2str(i)]);

            end

            
            filename= [sprintf('subj%02d_Regr_ER1', subi)]; % add one if all trials mixed by condition CHANGE TO Log
            save(filename,'B','R','stats','-v7.3');

%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% Compute Grandaverages of regressions here:

            cd(timeseries_folder);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(regression_folder);

            B_all= {};
            for subi=1:nSubjs;
                
                fname_Regr= sprintf('subj%02d_Regr_ER1',subi);
                pickupRegr(subi) = load(fname_Regr);
                
                templateER.avg = squeeze(pickupRegr(subi).B); %squeeze to collapse trials dimension
                B_all{subi}= templateER;
                
            end
            
            statsB_all= {};
            for subi=1:nSubjs;
                
                fname_Regr= sprintf('subj%02d_Regr_ER1',subi);
                pickupRegr(subi) = load(fname_Regr);
                
                templateER.avg = squeeze(pickupRegr(subi).stats); %squeeze to collapse trials dimension
                statsB_all{subi}= templateER;
                
            end
            

            cd(regression_folder);
            save('B_all_ER1', 'B_all', 'statsB_all','templateER','-v7.3');


            GAVG_B= ft_timelockgrandaverage([], B_all{1});
            
% Plot regressions:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.showlabels= 'yes';
            figure
            ft_multiplotER(cfg,GAVG_B);
            
%% 2b) Time-series: REGRESSION between RP amplitudes (Y) vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nTimes= 2001; % corresponds to [-3 +1]s

% 
% clear data_Path; %set right path
% data_Path= timeseries_folder; % need to fix it in the start up script
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in useful data:

    % Y= Time-locked amplitudes (alias RP avg)
    cd(timeseries_folder);
    load('pickupER');
    
%     for subi=1:nSubjs;
% 
%         fname_ER= sprintf('subj%02d_TimeS_bytrial',subi);
%         pickupER(subi) = load(fname_ER);
% 
%     end
    
    % X= Response times (alias Waiting Times)
    cd(behavioral_folder);
    load('pickupBehav');
    
%     for subi=1:nSubjs;
% 
%         fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
%         pickupBehav(subi) = load(fname_BehavData);
% 
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute regressions here:

cd(regression_folder);

        for subi=1:nSubjs;

            X= [pickupBehav(subi).RESPTIMES' ones(length(pickupBehav(subi).RESPTIMES'),1)];
            Y= pickupER(subi).avg_trl.trial; %.avg or .trial?
            if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

%             profile on

            clear B R stats;

            for i=1:nChans

                    for j= 1:nTimes

                       [B(:,i,j),~,R(:,i,j),~,stats(:,i,j)] = regress(Y(:,i,j),X);
            
                    end


                disp(['We are at channel' num2str(i)]);

            end


            filename= [sprintf('subj%02d_Regr_ER2', subi)]; % add one if all trials mixed by condition
            save(filename,'B','R','stats','-v7.3');
            
%             profile off;
%             profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of regressions here:

            cd(data_Path);
            
            load subj22_TimeS_one; % all conditions in one, trials averaged
            templateER= avg_one; % fake Fieldtrip structure for inserting correlation values
            
            % make a matrix with all the subjects
            cd(regression_folder);

            R_all= {};
            for subi=1:nSubjs;

                fname_Corr= sprintf('subj%02d_Regr_ER2',subi);
                pickupRegr(subi) = load(fname_Corr);

                templateER.avg = squeeze(pickupRegr(subi).R); %squeeze to collapse trials dimension
                B_all{subi}= templateER;

            end
            
            
            cd(regression_folder);
            save('B_all_ER2', 'B_all','templateER');


            GAVG_R= ft_timelockgrandaverage([], B_all{:});
            save 'GAVG_R_ER2' GAVG_R;
            
% Plot regressions:

            cfg=[];
            % cfg.preproc.lpfilter='yes';
            % cfg.preproc.lpfreq= 2; % Filter away alpha frequencies (8-12Hz): put 7 or 1 hz.
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 1.5;
            cfg.showlabels= 'yes';
            figure
            ft_multiplotER(cfg,GAVG_R);  

%% 3a) Time-frequency: correlation between power-spectra (Y) vs conditions (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nFreqs= 40;
nTimes= 61;
R= zeros(nChans,nFreqs,nTimes); P= zeros(nChans,nFreqs,nTimes);

data_Path= powerspectra_folder; % need to fix it in the start up script
cd(data_Path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Time-frequencies (alias TFR)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_TFR= sprintf('subj%02d_TFR_bytrial',subi); 
        pickupTFR(subi) = load(fname_TFR);
      
    end
    
%     save pickupTFRtrl pickupTFR;
    
    % X= Conditions (alias Time limits)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_Cond= sprintf('subj%02d_usefulinfo',subi);
        pickupCond(subi) = load(fname_Cond);
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlation here:

        
cd(regression_folder);
        
        for subi=1:nSubjs;

            profile on

             R(:)= 0; P(:)=0;

            pickupCond(subi).newcond(pickupCond(subi).newcond==Inf) = 32;  
            X= [pickupCond(subi).newcond' ones(length(pickupCond(subi).newcond'),1)];
            Y= pickupTFR(subi).TFR_trl.powspctrm;

                if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

                for i=1:nChans
                    for k= 1:nFreqs
                        for j= 1:nTimes

                            % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                            [tmpR,tmpP] = corrcoef(Y(:,i,k,j),X(:,1));
                            R(i,k,j)= tmpR(1,2);   P(i,k,j) = tmpP(1,2);

                        end
                    end

                    disp(['We are at channel' num2str(i)]);

                % save stuff
              
                filename= [sprintf('subj%02d_Corr_TFR1', subi)]; % add one if all trials mixed by condition
                save(filename,'R','P','-v7.3');

            end%     
            
            profile off;
            profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of correlations here:

            cd(data_Path);

            load subj22_TFR_oneCorrect; % all conditions in one, trials averaged
            templateTFR= TFR_one; % fake Fieldtrip structure for inserting correlation values

            % make a matrix with all the subjects
            cd(regression_folder);

            R_all= {};
            for subi=1:nSubjs;

                fname_Corr= sprintf('subj%02d_Corr_TFR1',subi);
                pickupCorrs(subi) = load(fname_Corr);

                templateTFR.powspctrm = squeeze(pickupCorrs(subi).R);
                R_all{subi}= templateTFR;


            end
            
            
            cd(regression_folder);
            save('R_all_TFR1', 'R_all','templateTFR');

            GAVG_R= ft_freqgrandaverage([], R_all{:});
            
% Plot correlations:

            cfg=[];
            % cfg.baseline= [-2 -1.5];
            % cfg.baselinetype= 'relative';
            cfg.zlim= [-0.05 -0.01];
%             cfg.xlim= [-2.3 -1.3];
            cfg.xlim= [-1.3 -0.7];
%             cfg.ylim= [11 31];
            cfg.ylim= [7 15];
            cfg.layout= 'eeg_64_NM20884N.lay';
            cfg.linewidth = 2;
            cfg.showlabels= 'yes';
            figure
%             ft_multiplotTFR(cfg,GAVG_R);
            ft_topoplotTFR(cfg,GAVG_R);
            
%% 3b) Time-frequency: correlation between power-spectra (Y) vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nFreqs= 40;
nTimes= 61;
R= zeros(nChans,nFreqs,nTimes); P= zeros(nChans,nFreqs,nTimes);

data_Path= powerspectra_folder; % need to fix it in the start up script
cd(data_Path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Time-frequencies (alias TFR)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_TFR= sprintf('subj%02d_TFR_bytrial',subi); 
        pickupTFR(subi) = load(fname_TFR);
        
    end
    
    % X= Response times (alias Waiting Times)
    cd(behavioral_folder);
    for subi=1:nSubjs;

        fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
        pickupBehav(subi) = load(fname_BehavData);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlation here:

cd(regression_folder);
        
        for subi=1:nSubjs;

            profile on

             R(:)= 0; P(:)=0;

          X= [pickupBehav(subi).RESPTIMES' ones(length(pickupBehav(subi).RESPTIMES'),1)];
            Y= pickupTFR(subi).TFR_trl.powspctrm;
            
                if isequal(length(Y(:,1)),length(X(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

                for i=1:nChans
                    for k= 1:nFreqs
                        for j= 1:nTimes

                            % [B(:,i,k,j),BINT,R] = regress(Y(:,i,k,j),X);
                            [tmpR,tmpP] = corrcoef(Y(:,i,k,j),X(:,1));
                            R(i,k,j)= tmpR(1,2);   P(i,k,j) = tmpP(1,2);
                            
                        end
                    end

                    disp(['We are at channel' num2str(i)]);

                % save stuff
              
                filename= [sprintf('subj%02d_Corr_TFR2', subi)]; % add one if all trials mixed by condition
                save(filename,'R','P','-v7.3');

                end
                
            profile off;
            profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of correlations here:

            cd(data_Path);

            load subj22_TFR_oneCorrect; % all conditions in one, trials averaged
            templateTFR= TFR_one; % fake Fieldtrip structure for inserting correlation values

            % make a matrix with all the subjects
            cd(regression_folder);

            R_all= {};
            for subi=1:nSubjs;

                fname_Corr= sprintf('subj%02d_Corr_TFR2',subi);
                pickupCorrs(subi) = load(fname_Corr);

                templateTFR.powspctrm = squeeze(pickupCorrs(subi).R);
                R_all{subi}= templateTFR;


            end
            
            
            cd(regression_folder);
            save('R_all_TFR2', 'R_all','templateTFR');
            

            GAVG_R= ft_freqgrandaverage([], R_all{:});
            
% Plot correlations:

            cfg=[];
            % cfg.baseline= [-2 -1.5];
            % cfg.baselinetype= 'relative';
            cfg.zlim= [-0.06 -0.02];
            cfg.xlim= [-2.75 -2.2];
            cfg.ylim= [12 31];
            % cfg.ylim= [15 30];
            cfg.layout = 'eeg_64_NM20884N.lay';
            cfg.linewidth = 2;
            cfg.showlabels= 'yes';
            figure
%             ft_multiplotTFR(cfg,GAVG_R);
            ft_topoplotTFR(cfg,GAVG_R);

            
%% 3c) Time-frequency: correlation between power-spectra (Y) sorted by conditions vs response times (X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters and right path:

nChans= 60;
nFreqs= 40;
nTimes= 61;
R= zeros(nChans,nFreqs,nTimes); P= zeros(nChans,nFreqs,nTimes);

data_Path= powerspectra_folder; % need to fix it in the start up script
cd(data_Path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in useful data:

    % Y= Time-frequencies (alias TFR)
    cd(powerspectra_folder);
    for subi=1:nSubjs;
        
        fname_TFR= sprintf('subj%02d_TFR_cond_trl',subi);
        pickupTFR(subi) = load(fname_TFR);
        
    end
    
    % X= Response times (alias Waiting Times)
    cd(behavioral_folder);
    for subi=1:nSubjs;

        fname_BehavData= sprintf('subj%02d_WaitingTimes',subi);
        pickupBehav(subi) = load(fname_BehavData);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute correlation here:

cd(regression_folder);

% sorted by condition
uconds=  [2     4     8    16   Inf];

        
        for subi=1:nSubjs;

            profile on

            R={}; P={};

            for condi = 1:length(uconds)

                X{condi}= [pickupBehav(subi).good_resps_cond{condi}' ones(length(pickupBehav(subi).good_resps_cond{condi}'),1)];
                Y{condi}= [pickupTFR(subi).TFR_cond{condi}.powspctrm];

                if isequal(length(Y{condi}(:,1)),length(X{condi}(:,1)))==1; disp('X & Y have correct dimensions'); else disp('X & Y DO NOT have compatible dimensions'); end;

                for i=1:nChans
                    for k= 1:nFreqs
                        for j= 1:nTimes

                            [tmpR,tmpP] = corrcoef(Y{condi}(:,i,k,j),X{condi}(:,1));
                            R{condi}(i,k,j) = tmpR(1,2);   P{condi}(i,k,j) = tmpP(1,2);

                        end
                    end

                    disp(['We are at channel' num2str(i)]);

                end

                % save stuff
                filename= [sprintf('subj%02d_Corr_TFR3', subi)]; % add one if all trials mixed by condition
                save(filename,'R','P','-v7.3');

                disp(['Condition ' num2str(condi) ' done']);

            end%

            profile off;
            profile viewer;
            disp(['Subject ' num2str(subi) ' done']);

        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sanity Check

for condi = 1:length(uconds)
    for i=1:nChans
        for k= 1:nFreqs
            for j= 1:nTimes
                if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at nTime n° ' num2str(j) ' of nFreq n°' num2str(j) ' nChan n° ' num2str(i) ' condi n° ' num2str(condi) ]); end;
            end
            if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at nFreq n° ' num2str(k) ' nChan n° ' num2str(i) ' condi n° ' num2str(condi)]); end;
        end
        if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at nChan n° ' num2str(i) ' condi n° ' num2str(condi)]); end;
    end
    if isequal(templateTFR{condi}.powspctrm(i,k,j),pickupCorrs(subi).R{condi}(i,k,j)); disp('YES'); else disp(['we are at condi n° ' num2str(condi)]); end;
end
        
% isequal(template{1}.powspctrm,template{2}.powspctrm,template{3}.powspctrm,template{4}.powspctrm,template{5}.powspctrm)
% isequal(pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi}, pickupCorrs(subi).R{condi})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
% Compute Grandaverages of correlations here:

            cd(data_Path);

            load subj22_TFR_cond; % all conditions in one, trials averaged
            templateTFR_cond= TFR_cond; % fake Fieldtrip structure for inserting correlation values

            % make a matrix with all the subjects
            cd(regression_folder);

           
            R_all= {};
            for subi=1:nSubjs;
              
                fname_Corr= sprintf('subj%02d_Corr_TFR3',subi);
                pickupCorrs(subi) = load(fname_Corr);
                
                for condi = 1:length(uconds)
                    
                    templateTFR_cond{condi}.powspctrm= pickupCorrs(subi).R{condi}; 
                    R_all{subi,condi}= templateTFR_cond{condi};
                    
                end
                
            end

            save 'R_all_TFR3' 'R_all';

%Grandaverage

            GAVG_R= {};
            for condi = 1:length(uconds)

                GAVG_R{condi}= ft_freqgrandaverage([], R_all{:,condi});

            end
            
            save GAVG_R_TFR3 GAVG_R;

% Normalize grandavg values

            for condi = 1:length(uconds)

                GAVG_R{condi}.powspctrm = normalize(GAVG_R{condi}.powspctrm);

            end
            
            save GAVG_R_TFR3_norm GAVG_R;


% Plot correlations:

            for condi = 1:length(uconds)

                cfg=[];

                % cfg.baseline= [-2 -1.5];
                % cfg.baselinetype= 'relative';
                % cfg.zlim= [-3e-25 3e-25];
                cfg.xlim= [-2.75 0];
                cfg.ylim= [0 30];
                cfg.zlim= [-0.2 0.2];
                % cfg.ylim= [15 30];
                cfg.layout = 'eeg_64_NM20884N.lay';
                cfg.linewidth = 2;
                cfg.showlabels= 'yes';
                figure

                ft_multiplotTFR(cfg,GAVG_R{condi});

            end

