%% Mean/Median RP amplitude for a specific time window
% ======================================================================= %
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on July 2018.
% EXPERIMENT: Timelimit_2018


%{ 

    SCOPE: Script for extracting mean RP amplitude for plotting and stats.

    OUTPUT: datamatrix_premov_20{i,k}, mean_premov_amp_ch20(i,k),
          median_premov_amp_ch20(i, k).
    
    HOW:
    
        1)Load Timelimit_*_subj**_EEG_clean_concat_rej_interp.mat.
        2)Create a new matrix for subjects (11) x conditions (5) and extraxt amp.
        3)Do mean and median of datamatrix_premov_**{i,k}.avg.
        4)Save datamatrix_premov_**, mean_premov_amp_ch**,median_premov_amp_ch**.
        5)Plot mean_all, median_all.

%}

% FIX ME: the loop to load each preprocessed file still not working

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

%% Consider only GOOD SUBJECTS 

% With clear Readiness Potential(= negative slope vs flat slope or positive slope).
good_subjs = [2 3 5 6 7 8 10 11 13 15 17 18 20]; % % removed: subj04,subj14,subj19; remove subj 12 & 16 for channel 28 and 30.
    nGoodSubjs = length(good_subjs);

%% Load timelocked averages from each subject's folders

for i=1:nSubjs;
    
    if i==1
%         fname = sprintf('subj%02d_noHPI/Timeseries/Amplitudes/Baseline_Haggard/timeseries_EEG.mat',good_subjs(i)); %nGoodSubjs(i)
        fname = sprintf('subj%02d_noHPI/Timeseries/Amplitudes/No_Baseline/timeseries0_EEG.mat',i);
    else
%         fname = sprintf('subj%02d/Timeseries/Amplitudes/Baseline_Haggard/timeseries_EEG.mat',i);
        fname = sprintf('subj%02d/Timeseries/Amplitudes/No_Baseline/timeseries0_EEG.mat',i);
    end
    if isfile(fname)
        pickupSub(i) = load(fname);
    end
end

%% Create a new matrix for subjects (11) x conditions (5) and extraxt amp
% FIX ME IN A SMART WAY
%=========================================================================%
% channel 20 (29 Oct: added 2 subjects!!)
% datamatrix_premov_20 = [];
% mean_premov_amp_ch20 = [];
% median_premov_amp_ch20= [];

% % fix this
% good_subjects = [2 3 5 6 7 8 10 11 12 13 15 16 17]; % removed: subj04,subj09,subj14.
% nGoodSubjects = length(good_subjects);
datamatrix_premov= struct('ch20',[],'ch28', [],'ch30',[],'ROI', []); 
mean_premov_amp= struct('ch20',[],'ch28', [],'ch30',[],'ROI', []); 
% sem_premov_amp= struct('ch20',[],'ch28', [],'ch30',[]);

for i= 1:nSubjs; %nGoodSubjects or nSubjects
    
    for k= 1:5
        
% %         avg2matrix{i,k}= pickupSub(i).avg{k}; %avg_EEG
% %         stdmatrix{i,k}= pickupStd(i).across_stdev{k}.avg
%         cfg = [];
%         cfg.latency = [-1 -.2];
%         cfg.channel = 'EEG020';
%         datamatrix_premov.ch20{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
%         mean_premov_amp.ch20(i, k) = mean(datamatrix_premov.ch20{i,k}.avg);
% %         sem_premov_amp.ch20(i, k) = sem(mean_premov_amp.ch20(i, k),1);
%         
%         cfg = [];
%         cfg.latency = [-1 -.2];
%         cfg.channel = 'EEG028';
%         datamatrix_premov.ch28{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
%         mean_premov_amp.ch28(i, k) = mean(datamatrix_premov.ch28{i,k}.avg);
% %         sem_premov_amp.ch28(i, k) = sem(mean_premov_amp.ch28(i, k),1);
% 
        cfg = [];
        cfg.latency = [-1 0]; %[-1 -.2]
        cfg.channel = 'EEG030';
        datamatrix_premov.ch30{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
        mean_premov_amp.ch30(i, k) = mean(datamatrix_premov.ch30{i,k}.avg);
%         sem_premov_amp.ch30(i, k) = sem(mean_premov_amp.ch30(i, k),1);

%         cfg = [];
%         cfg.latency = [-1 0];%[-1 -.2]
%         cfg.channel = {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'};
%         datamatrix_premov.ROI{i,k} = ft_selectdata(cfg, avgmatrix{i,k});
%         datamatrix_premov.ROI{i,k}= mean(datamatrix_premov.ROI{i,k}.avg,1);
%         mean_premov_amp.ROI(i, k) = mean(datamatrix_premov.ROI{i,k});
        
    end
end

%% Save

% Create the folder if it doesn't exist already.
if input('Save MEAN AMPLITUDE results? RISK OF OVERWRITING  (1/0) ... ')
    
    results_folder= [parent_folder, '/Results'];
         if ~exist(fullfile(parent_folder)); mkdir(fullfile(results_folder)); end;
    cd(results_folder);

   save datamatrix_premov datamatrix_premov; save mean_premov_amp mean_premov_amp 
   
end

%% PLOT (24 Oct 2018)
cd(resfolder);

load('mean_premov_amp')


figfolder = fullfile(parent_folder,'/Figures');
cd(figfolder);
%%
%Function rule for Recursive sequence (24 Oct)
a=2; r=2;n=5;
s = a*r.^(0:n-1); 

mean_all= struct('ch20',[],'ch28', [],'ch30',[],'ROI', []);
sem_all= struct('ch20',[],'ch28', [],'ch30',[]),'ROI', [];

% for k=1:5
    
    mean_all.ch20= mean(mean_premov_amp.ch20(:,:));
    sem_all.ch20= sem(mean_premov_amp.ch20(:,:),1);
    
    mean_all.ch28= mean(mean_premov_amp.ch28(:,:));
    sem_all.ch28= sem(mean_premov_amp.ch28(:,:),1);
    
    
    mean_all.ch30= mean(mean_premov_amp.ch30(OKsubjs,:));
    sem_all.ch30= sem(mean_premov_amp.ch30(OKsubjs,:),1);
    
    mean_all.ROI= mean(mean_premov_amp.ROI(OKsubjs,:));
    sem_all.ROI= sem(mean_premov_amp.ROI(OKsubjs,:),1);
    
% end

% f2=figure(1)
figure
errorbar(s,mean_all.ch20,sem_all.ch20,'-or','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean')  % RED= mean
hold on
errorbar(s,mean_all.ch28,sem_all.ch28,'-ob','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean')  % RED= mean
hold on
errorbar(s,mean_all.ch30,sem_all.ch30,'-og','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean')  % RED= mean
set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Conditions (sec)')
ylabel('Mean RP amplitudes (\muV')
title(['Mean RP voltage for ' num2str(i) ' subjs, all channels'])
legend('Channel 20','Channel 28', 'Channel 30','show','Location','best')
xlim([0 Inf])
hold on

% filename= ['RP_'];
filename1= ['RPsLog_amps_chan30_all.png'];
saveas(f1,filename1)
f2=figure(2);
filename2= ['Residual_fit_EEG_chan30_all.png'];
saveas(f2,filename2)

%% barplot
figure
bar(s,mean_all.ch20)
hold on
errorbar(s,mean_all.ch20,sem_all.ch20,'r.','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean')
set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Conditions (sec)')
ylabel('RP amplitudes (\muV')
title('Channel FC1')

hold off
bar(s,mean_all.ch28)
hold on
errorbar(s,mean_all.ch28,sem_all.ch28,'r.','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean')
set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Conditions (sec)')
ylabel('RP amplitudes (\muV')
title('Channel C3')

hold off
bar(s,mean_all.ch30)
hold on
errorbar(s,mean_all.ch30,sem_all.ch30,'r.','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean')
set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Conditions (sec)')
ylabel('RP amplitudes (\muV')
title('Channel Cz')

bar(s,mean_all.ROI)
hold on
errorbar(s,mean_all.ROI,sem_all.ROI,'r.','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean')
set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Conditions (sec)')
ylabel('RP amplitudes (\muV')
title('ROI: channels 20,21,29,30,31,39,40')

%% prepare matrix for stats
% let's remove the 3rd condition but handle with care (if you run it twice
% it will remove 2 columns!)
%mean_premov_amp(:,3) = [];

% column 1: dependent measure
% Xamps= mean_premov_amp(:);
% 
% %Column2: IV, Independent variable (conditions, 1:5 x17)
% Xcond= repmat([1:nConds],nGoodSubjects,1)
% 
% %Column3: Subjects, 11
% Xsubj = repmat([1:nGoodSubjects]',1,nConds)
% 
% % MATRIX 
% X = [Xamps(:) Xcond(:) Xsubj(:)]
% 
% %% now run anova
% 
% [x,y,z]=RMAOV1(X)
% 
% % F(4,10) = 2.951, p=0.0315
% 
% %% extract values for post hoc t-tests
% 
% [H54,P54,CI54,STATS54]= ttest(mean_premov_amp(:,5), mean_premov_amp(:,4))
% 
% [H43,P43,CI43,STATS43]= ttest(mean_premov_amp(:,4), mean_premov_amp(:,3))
% 
% [H32,P32,CI32,STATS32]= ttest(mean_premov_amp(:,3), mean_premov_amp(:,2))
% 
% [H21,P21,CI21,STATS21]= ttest(mean_premov_amp(:,2), mean_premov_amp(:,1))
% 
% 
% [H53,P53,CI53,STATS53]= ttest(mean_premov_amp(:,5), mean_premov_amp(:,3))
% 
% [H52,P52,CI52,STATS52]= ttest(mean_premov_amp(:,5), mean_premov_amp(:,2))
% 
% [H51,P51,CI51,STATS51]= ttest(mean_premov_amp(:,5), mean_premov_amp(:,1))
% 
