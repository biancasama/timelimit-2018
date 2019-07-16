%% Plot time-series between (within) subject
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: re_created on June-July 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: 

    OUTPUT: 

    FIXME: **

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
    
%% load stuff 

cd(timeseries_folder);
load Grand_ER_Ind;

%% EELKE's way (http://eelkespaak.nl/blog/customizing-common-m-eeg-plots-part-1-the-event-related-potential-field-erp-f/)

figure('color','white');
for condi=1:5
    
    tl = Grand_ER{condi};
    % Find the index of the channels of interest in the data
    chaninds = match_str(tl.label, {'EEG030'});
    
    % Compute average over channels
    erp = tl.avg(chaninds,:); % add squeeze if indiv subjs | tl.individual
    
    % Plot
    %     figure();
    hold on;
    plot(tl.time, erp);
    
    micro = 1e-6;
    erp_microV = erp ./ micro;
    
    %     figure();
    %     colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]};
    Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]};
    plot(tl.time, erp_microV,'LineWidth',2.5); % 'Color',Colors{condi}
    %     xlim([-3.0 1.0]);
    xlabel('Time (s)');
    ylabel('Amplitude electric potential (\mu V)');
    
    %grid
    
    % Change axis
    %     ax = gca(); % this Gets the Current Axis so we can set properties
    %     ax.XAxisLocation = 'origin';
    %     ax.YAxisLocation = 'origin';
    %     %ax.TickDir = 'out';
    %
    %     % Remove the box around the plot, while we're at it:
    %     box off;
    %
    %     % And move the x-axis label to underneath the axis:
    %     ax.XLabel.Position(2) = -30;
    clear tl erp
end
%% EELKE's way + boundedline = VERSION USED 

addpath(genpath('/Users/bt_neurospin/matlab/kakearney-boundedline-pkg-50f7e4b'));
conds= [2 4 8 16 Inf];
mycolormap= [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.4940 0.1840 0.5560]; %linespec
channels= [20 28 30]; %20= FC1; 28= C3; 30= CZ (EEG cap 60 electrodes)
chan_names= {'FC1','C3','Cz'};
ROI = [20 21 29 30 31 39 40];
chaninds= [30 20 28 {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'}];
if chaninds == 30
    chan_names= {'Cz'}
elseif chaninds == 20
    chan_names= {'FC1'}
elseif chaninds == 28
    chan_names== {'C3'}
else chaninds== {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'}; %Undefined operator '==' for input arguments of type 'cell'.
    chan_names= {'ROI'}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f1 is for ALL SUBJS
    f1=figure('color','white');
    for condi=1:length(conds);
        %     figure('color','white');
        tl = Grand_ER{condi};
%         colors= Colors{condi};

        % Find the index of the channels of interest in the data
        %chaninds = match_str(tl.label, {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'});
        chaninds = match_str(tl.label, {'EEG030'});
        % Compute average over channels
        erp = squeeze(mean(tl.individual(:,chaninds,:),1)); % add squeeze if indiv subjs | tl.individual
        micro = 1e-6;
        erp_microV = erp ./ micro;
        
        % First, get the individual trial data, averaged over our channels of
        % interest:
        trialdata = squeeze(mean(tl.individual(:,chaninds,:), 2)) ./ micro;
        % Use the standard deviation over trials as error bounds:
        bounds = sem(trialdata, 1); %sem or std?
        % boundedline will replace the call to plot():
        b=boundedline(tl.time, erp_microV, bounds, 'cmap', mycolormap(condi,:),'alpha'); % alpha makes bounds transparent
        b.LineWidth= 2;
        
        % If you want to add effects, do it here:
        %stattimes= [-1.0320 -0.316];
%         stattimes1= linspace(-0.742,-0.26,60);
% %         stattimes2= linspace(-1.398,-1.098,40);
%         y1= linspace(1,1,60);
% %         y2= linspace(1,1,40);
%         txt1 = {'p < .05'};
%         t1=text(-0.75,2,txt1,'FontSize',25);
% %         txt2 = {'p < .05'};
% %         t2=text(-1.4,2,txt2,'FontSize',25);
%         l1=line(stattimes1,y1); %,'Color','black','Marker','*'
%         l1.Color= 'black'; 
%         l1.LineWidth= 3;
%         l1.Marker= '*';
%         l1.MarkerSize= 10;
%         darkGrey1  = [0.4 0.4 0.4]; darkGrey2  = [0.2 0.2 0.2]; darkGrey3  = [0.1 0.1 0.1]; 
% %         l2=line(stattimes2,y2); %,'Color','black','Marker','*'
% %         l2.Color= darkGrey1;
% %         l2.LineWidth= 2;
% %         l2.Marker= '*';
% %         l2.MarkerSize= 10;
        
        xlabel('Time (s)','FontSize',34);
        ylabel('ERP amplitude (\muV)','FontSize',34);
        title(cellstr(chan_names),'FontSize',34);
        set(gca,'FontSize',34);
        %     legend(subset,'FontSize',12)
        %     ax = gca();
        %     ax.XAxisLocation = 'origin';
        %     ax.YAxisLocation = 'origin';
        %     ax.TickDir = 'out';
        %     box off;
        %     ax.XLabel.Position(2) = -60;
        hold on;
        hold off;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f2 is only for OKSubjs
    f2=figure('color','white');
    for condi=1:length(conds);
        %     figure('color','white');
        tl = Grand_ER{condi};

        % Find the index of the channels of interest in the data
        chaninds = match_str(tl.label, {'EEG030'});
        
        % Compute average over channels
        erp = squeeze(mean(tl.individual(:,chaninds,:),1)); % add squeeze if indiv subjs | tl.individual
        micro = 1e-6;
        erp_microV = erp ./ micro;
        
        % First, get the individual trial data, averaged over our channels of
        % interest:
        trialdata = squeeze(mean(tl.individual(:,chaninds,:), 2)) ./ micro;
        
        % Use the standard deviation over trials as error bounds:
        bounds = sem(trialdata, 1); %sem or std?
        % boundedline will replace the call to plot():
        b=boundedline(tl.time, erp_microV, bounds, 'cmap', mycolormap(condi,:),'alpha'); % alpha makes bounds transparent
        b.LineWidth= 2;
        
        % If you want to add effects, do it here:
        stattimes1= linspace(-0.742,-0.26,60);%linspace(-0.722,-0.448,40);
        %stattimes2= linspace(-1.868,-1.63,20);
        y1= linspace(1,1,60);
%         y2= linspace(1,1,20);
        txt1 = {'p < .05'};
        t1=text(-0.75,2,txt1,'FontSize',25);
%         txt2 = {'p < .05'};
%         t2=text(-1.87,4,txt2,'FontSize',25);
        l1=line(stattimes1,y1); %,'Color','black','Marker','*'
        l1.Color= 'black'; 
        l1.LineWidth= 3;
        l1.Marker= '*';
        l1.MarkerSize= 10;
        darkGrey1  = [0.4 0.4 0.4]; darkGrey2  = [0.2 0.2 0.2]; darkGrey3  = [0.1 0.1 0.1]; 
%         l2=line(stattimes2,y2); %,'Color','black','Marker','*'
%         l2.Color= darkGrey2;
%         l2.LineWidth= 2;
%         l2.Marker= '*';
%         l2.MarkerSize= 10;
        
        xlabel('Time (s)','FontSize',34);
        ylabel('ERP amplitude (\muV)','FontSize',34);
        title(cellstr(chan_names),'FontSize',34);
        set(gca,'FontSize',34);
%         lgd= legend(,'FontSize',12)
        %     ax = gca();
        %     ax.XAxisLocation = 'origin';
        %     ax.YAxisLocation = 'origin';
        %     ax.TickDir = 'out';
        %     box off;
        %     ax.XLabel.Position(2) = -60;
        hold on;
        hold off;
    end
    
%% per individual subjects (code not working)
nSubjs=22;
condi=5;
    for subi=1:length(nSubjs);
        figure('color','white');
        tl = Grand_ER{condi};

        % Find the index of the channels of interest in the data
        %chaninds = match_str(tl.label, {'EEG020','EEG021','EEG029','EEG030','EEG031','EEG039','EEG040'});
        chaninds = match_str(tl.label, {'EEG030'});
        
        % Compute average over channels
        erp = squeeze(tl.individual(subi,chaninds,:)); % add squeeze if indiv subjs | tl.individual
        micro = 1e-6;
        erp_microV = erp ./ micro;
        
        % First, get the individual trial data, averaged over our channels of
        % interest:
        trialdata = (tl.individual(subi,:) ./ micro);
        
        % Use the standard deviation over trials as error bounds:
        bounds = std(trialdata, 1); %sem or std?
        % boundedline will replace the call to plot():
        b=boundedline(tl.time, erp_microV, bounds, 'cmap', mycolormap(condi,:),'alpha'); % alpha makes bounds transparent
        b.LineWidth= 2;
        
        % Add all our previous improvements:
        xlabel('Time (s)','FontSize',34);
        ylabel('ERP amplitude (\muV)','FontSize',34);
        title(['subj ' num2str(subi) cellstr(chan_names)],'FontSize',34);
        set(gca,'FontSize',34);
        %     legend(subset,'FontSize',12)
        %     ax = gca();
        %     ax.XAxisLocation = 'origin';
        %     ax.YAxisLocation = 'origin';
        %     ax.TickDir = 'out';
        %     box off;
        %     ax.XLabel.Position(2) = -60;
%         hold on;
        %filename= ['RPprofile_subj_' num2str(subi) '.png'];
%         hold off;
    end
    
%% OLD way:plot to check if it shows RP

conditions_all= [1 2 3 4 5];
n_conditions= length(conditions_all);
channels= [20 28 30]; %20= FC1; 28= C3; 30= CZ (EEG cap 60 electrodes)

for j= 1:numel(channels)
    h(j)=figure;
    for condi= 5
        
        selected_data = Grand_Avg{condi}.avg(channels(j),:); % 1001:1500 samples for statistics
        selected_time = Grand_Avg{condi}.time(:);
        h(j,condi)=figure;
        % subplots?
        %         subplot(3,5,15)
        plot(selected_time, selected_data,'LineWidth',3)
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        title(['Grandaverage Readiness potential from condition ' num2str(condi) ', channel ' num2str(channels(j))])
        %         title(['Cond ' num2str(condi) ', chan' num2str(channels(j))])
        xlim([-3.0 1.0])
        %         ylim([-1e-13 3e-13]) %??
        %         hold on
        % filename= ['Grandavg_RP_allconds', '_Chan_' int2str(channels(j)) '.png'];
        
    end
    filename= ['Grandavg_RP_Inf', '_Chan_' int2str(channels(j)) '.png'];
%         saveas(h(j),filename); % you will save 15 figures (5x3)
        hold on
        hold off
end

% topoplot (try for Early RP [-1 -0.2] vs Late RP [-0.2 0])
%% OLD way: regular plot

conditions_all= [1 2 3 4 5];
n_conditions= length(conditions_all);
channels= [20 28 30]; %20= FC1; 28= C3; 30= CZ (EEG cap 60 electrodes)
chan_names= {'FC1','C3', 'Cz'};

for j= 1:numel(channels)
    k(j)=figure;
    for condi= 1: n_conditions;
        
        colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]} % Cell array of colors.
        hold on
        p1=plot(Grand_Avg{condi}.time(:),Grand_Avg{condi}.avg(channels(j),:),'LineWidth',2); %this time use samples on the time axes
        hold on
        title(['Grandaverage RP amplitude, N= ' int2str(numel(OKsubjs)), ', Channel ' chan_names{j} ]);
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        legend('2sec','4sec','8ec','16sec','Inf','Location','SouthEast')% COLORS DO NOT MATCH
        hold on
        filename= ['RP_grandavg_Chan_' int2str(channels(j)) '_allconds', '.png'];
        
    end
%     saveas(k(j),filename); %HOW TO SAVE ONLY LAST FIGURE FROM A LOOP??
    hold on
    hold off
end

channels= [20,21,29,30,31,39,40];
chan_names= {'FC1','C3', 'Cz'};

figure;
    for condi= 1: n_conditions;
        
        colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]} % Cell array of colors.
        hold on
        p1=plot(Grand_Avg{condi}.time(:),Grand_Avg_ROI{condi},'LineWidth',2); %this time use samples on the time axes
        hold on
        title(['Grandaverage RP amplitude, N= ' int2str(numel(OKsubjs)), ', ROI channels' ]);
        xlabel('Time (s)');
        ylabel('Amplitude (\muV)');
        legend('2sec','4sec','8ec','16sec','Inf','Location','SouthEast')% COLORS DO NOT MATCH
        hold on
        filename= ['RP_grandavg_ROI_chans.png'];
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLot mean amplitude averaged for specific window of time
% load stuff created with extract_mean-amp code

cd(timeseries_folder);
load 'mean_all' 'sem_all';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END