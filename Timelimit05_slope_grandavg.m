%% SLOPE averages for Timelimit2018
%=========================================================================%
% This routine computes the MEAN slopes' regressions across
% participants and plots them vs amplitudes' grandaverages computed through
% timelockgrandaverage.

% AUTHOR: Bianca Trov?.
% DATE: 25 July 2018; modified: 25 September 2018.

% Output: Slope_all, only mean slopes values, saved structure looks like
%                   Slope_all.ch20{condi}(subj)
%                   Slope_all.ch28{condi}(subj)
%                   Slope_all.ch30{condi}(subj)
%         Intercept_all, only mean intercept values, saved structure looks like
%                   Intercept_all.ch20{condi}(subj)
%                   Intercept_all.ch28{condi}(subj)
%                   Intercept_all.ch30{condi}(subj)
%         Y_all, mean regression line, saved structure looks like
%                   Y_all.ch20{condi}(subj,:)
%                   Y_all.ch28{condi}(subj,:)
%                   Y_all.ch30{condi}(subj,:)
%         GrandSlopeallmatrix, redundant, containing
%                   = [GrandSlopeall.ch20;GrandSlopeall.ch28;GrandSlopeall.ch30]
%         GrandYallmatrix, for plotting purposes, containing
%                   = [GrandYall.ch20;GrandYall.ch28;GrandYall.ch30]

% HOW:
%{ 

Load all the files called SlopeRP.mat (1=Slope+2=Intercept)
Take the mean for each channel
SlopeRP{condi}(1,channel)
Save GrandSlopeRP
Do the same for Y
Save GrandY
Plot it against grandaveraged RP? Or plot only the slopes?

%}

% FIX ME: Plots part: legend.
% FIX ME: Compute error??
%=========================================================================%
%% Set path (27.09)
clearvars;
% Define some paths 
if strcmp(computer, 'MACI64')% on my laptop
    script_Path= '/Volumes/USB_DISK/TIMELIMIT_backup/SCRIPTS_ANALYSES/MEEG'; % here you find all the scripts for preprocessing/analysing MEEG data
    data_Path = '/Volumes/USB_DISK/TIMELIMIT_backup/MEEG_fif_files'; % = parent_folder: all the raw data are stored here.
    ft_Path = '//Users/bt_neurospin/matlab/FIELDTRIP/fieldtrip-20170405'; % Fieldtrip tools 
    tool_Path= '/Users/bt_neurospin/matlab'; % other useful functions for matlab (written by Aaron)
end

% add some paths
addpath(genpath(script_Path)); % general pre-proc path
addpath(genpath(tool_Path));
addpath(ft_Path); ft_defaults; % start up fieldtrip [NEW]
addpath([ft_Path '/engine']); % start up FT engines [NEW] NOTE: what is this actually doing? 

parent_folder=data_Path; 
subj_folders = dir(fullfile(parent_folder, 'subj*'))
nSubjects= length(subj_folders);
conditions_all= [1 2 3 4 5];
n_conditions= length(conditions_all);
channels= [20 28 30];
good_subjects = [2 3 5 6 7 8 10 11 13 15 17]; % % removed: subj04,subj14,subj19; remove subj 12 & 16 for channel 28 and 30.
nGoodSubjects = length(good_subjects);

%% Load data: SLOPE
% CORRECTION here (25.09.2018): Intercept_all > Slope_all
% 13th of November: added the other participants that weren't good too.

Slope_all= struct('ch20',[],'ch28', [],'ch30',[]); % Best_Slope_all
for subj= 1: nGoodSubjects % nGoodSubjects
    for condi= 1: n_conditions
        current_subject = fullfile(parent_folder, subj_folders(subj).name) % subj_folders(good_subjects(subj)).name)
        load([current_subject,'/Timeseries/Slopes/Baseline_Haggard/SlopeRP.mat']) %'/Slopes/SlopeRP.mat'
        Slope_all.ch20{condi}(subj)= SlopeRP{condi}(1,1); % 1= ch 20
        Slope_all.ch28{condi}(subj)= SlopeRP{condi}(1,2); % 2= ch 28
        Slope_all.ch30{condi}(subj)= SlopeRP{condi}(1,3); % 3= ch 30
       
    end
end

%% Load data: INTERCEPT
% Take 2 from SlopeRP for Intercept

% Intercept_all= struct('ch20',[],'ch28', [],'ch30',[]);
% for subj= 1:nGoodSubjects
%     for condi= 1: n_conditions
%         current_subject = fullfile(parent_folder, subj_folders(good_subjects(subj)).name)
%         load([current_subject,'/Timeseries/Slopes/Baseline_Haggard/Slopes/SlopeRP.mat'])
%         Intercept_all.ch20{condi}(subj)= SlopeRP{condi}(2,1); % 1= ch 20
%         Intercept_all.ch28{condi}(subj)= SlopeRP{condi}(2,2); % 2= ch 28
%         Intercept_all.ch30{condi}(subj)= SlopeRP{condi}(2,3); % 3= ch 30
%        
%     end
% end

 Intercept_all= struct('ch20',[],'ch28', [],'ch30',[]);
for subj= 1:nGoodSubjects
    for condi= 1: n_conditions
        current_subject = fullfile(parent_folder, subj_folders(subj).name)
        load([current_subject,'/Timeseries/Slopes/Baseline_Haggard/SlopeRP.mat']) %'/Slopes/SlopeRP.mat'
        Intercept_all.ch20{condi}(subj)= SlopeRP{condi}(2,1); % 1= ch 20
        Intercept_all.ch28{condi}(subj)= SlopeRP{condi}(2,2); % 2= ch 28
        Intercept_all.ch30{condi}(subj)= SlopeRP{condi}(2,3); % 3= ch 30
       
    end
end

%% Load data: REGRESSION LINE

Y_all= struct('ch20',[],'ch28', [],'ch30',[]);
for subj= 1:nGoodSubjects
    for condi= 1: n_conditions
        current_subject = fullfile(parent_folder, subj_folders(subj).name)
        load([current_subject,'/Timeseries/Slopes/Baseline_Haggard/Y.mat'])
        Y_all.ch20{condi}(subj,:)= Y{condi}(1,:);
        Y_all.ch28{condi}(subj,:)= Y{condi}(2,:);
        Y_all.ch30{condi}(subj,:)= Y{condi}(3,:);
       
    end
end

%% Compute mean & SEM
% stderror = std( data ) / sqrt( length( data ))

for condi= 1: n_conditions
    
    GrandSlopeall.ch20{condi}=mean(Slope_all.ch20{condi}(:))
%     SEMSlopeall.ch20{condi}=std((Slope_all.ch20{condi}(:)),sqrt(length(Slope_all.ch20{condi}(:))))
    GrandSlopeall.ch28{condi}=mean(Slope_all.ch28{condi}(:))
    GrandSlopeall.ch30{condi}=mean(Slope_all.ch30{condi}(:))
    
end

GrandSlopeallmatrix= [GrandSlopeall.ch20;GrandSlopeall.ch28;GrandSlopeall.ch30]

%% Y (regression line)
for condi= 1: n_conditions
    
    GrandYall.ch20{condi}=mean(Y_all.ch20{condi}(:,1:500))
%     SEMYall.ch20{condi}=std((Y_all.ch20{condi}(:)),sqrt(length(Y_all.ch20{condi}(:))))
    GrandYall.ch28{condi}=mean(Y_all.ch28{condi}(:,1:500))
    GrandYall.ch30{condi}=mean(Y_all.ch30{condi}(:,1:500))
    
end

GrandYallmatrix= [GrandYall.ch20;GrandYall.ch28;GrandYall.ch30]

%% Save stuff in the right folder (25.09.18)

% Create the folder if it doesn't exist already.
% [status, msg] = mkdir('Results'); %IT'S NOT WORKING PROPERLY
% 
% cd(fullfile(parent_folder,'/Results'));

%create folder if it doesn't already exist
[status, msg] = mkdir(parent_folder,'Results/Timeseries/Slopes');
cd(fullfile(parent_folder,'Results/Timeseries/Slopes'));

save Slope_all Slope_all;
save Intercept_all Intercept_all;
save Y_all Y_all;
save GrandSlopeallmatrix GrandSlopeallmatrix;
save GrandYallmatrix GrandYallmatrix;

%% Plot all the conditions in a unique plot
% Modif. added 25.09.18
load Grand_Avg % from (parent_folder,'/Results')
figfolder = fullfile(parent_folder,'/Figures');
cd(figfolder);
[status, msg] = mkdir('Finalfig');
cd(fullfile(pwd,'/Finalfig'));
[status, msg] = mkdir('Slopes_fig');
cd(fullfile(pwd,'/Slopes_fig'));


for j= 1:numel(channels)
    h(j)=figure;
    for condi= 1: n_conditions;
        
        p=plot(Grand_Avg{condi}.time(:),Grand_Avg{condi}.avg(channels(j),:),'LineWidth',1) %this time use samples on the time axes
%         legend('2 sec','4 sec','8sec','16sec','Inf','Location','Northwest','off')
        hold on
        s=plot([-1:(1/500):-0.0020],GrandYallmatrix{j,condi},'LineWidth',2.5)
        legend('Slope 2 sec','Slope 4 sec','Slope 8 sec','Slope 16 sec','Slope Inf','Location','Best')
        hold on
        
        title(['GrandReadiness Potential for all conditions, channel ' int2str(channels(j)) ])
        xlabel('Time (s)')
        ylabel('Amplitude (\muV)')
        
        % LEGEND NEEDED
    end
    filename= ['GrandSlopesRP_allconds', '_Chan_' int2str(channels(j)) '.png'];
    saveas(h(j),filename); % you will save 15 figures (5x3)
    hold on
    hold off
end

%% Plot only channel 20
% LEGEND NOT WORKING - what about removing useless loop for channels?
j=1;
k1=figure;
hold on
for condi= 1: n_conditions;
    
    colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]} % Cell array of colors.
    hold on
    p1=plot(Grand_Avg{condi}.time(:),Grand_Avg{condi}.avg(channels(j),:)); %this time use samples on the time axes
    hold on
    s1(condi)=plot([-1:(1/500):-0.0020],GrandYallmatrix{j,condi},'LineWidth',2,'Color',colors{condi});
    title(['GrandRP Slope fit for all conditions, ' ' Channel ' int2str(channels(j)) ]);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
%     legend('2sec','4sec','8ec','16sec','Inf','Location','Best')% COLORS DO NOT MATCH
    hold on
    filename= ['RP_grandslope_Chan_' int2str(channels(j)) '_cond_' int2str(condi), '.png'];
    saveas(k1,filename); %HOW TO SAVE ONLY LAST FIGURE FROM A LOOP??
    hold on
    hold off
end

%  Legend=cell(2,1)%  two positions 
%  Legend{1}=' Your data 1' ;
%  Legend{2}=' Your data 2';
%  legend(Legend);
%% Plot only channel 28
j=2;
k2=figure;
hold on
for condi= 1: n_conditions;
    
    colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]} % Cell array of colors.
    hold on
    p2=plot(Grand_Avg{condi}.time(:),Grand_Avg{condi}.avg(channels(j),:)); %this time use samples on the time axes
    hold on
    s2(condi)=plot([-1:(1/500):-0.0020],GrandYallmatrix{j,condi},'LineWidth',2,'Color',colors{condi});
    title(['GrandRP Slope fit for all conditions, ' ' Channel ' int2str(channels(j)) ]);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    %legend([s1],'2sec','4sec','8ec','16sec','Inf','Location','Best')% COLORS DO NOT MATCH
    hold on
    filename= ['RP_grandslope_Chan_' int2str(channels(j)) '_cond_' int2str(condi), '.png'];
    saveas(k2,filename); %HOW TO SAVE ONLY LAST FIGURE FROM A LOOP??
    hold on
    hold off
end

% LEGEND NOT WORKING

%% Plot only channel 30

j=3;
k3=figure;
hold on
for condi= 1: n_conditions;
    
    colors = {[0 0 1],[1 0 0],[0 1 0],[0 0 0],[1 1 0]} % Cell array of colors.
    hold on
    p3=plot(Grand_Avg{condi}.time(:),Grand_Avg{condi}.avg(channels(j),:)); %this time use samples on the time axes
    hold on
    s3(condi)=plot([-1:(1/500):-0.0020],GrandYallmatrix{j,condi},'LineWidth',2,'Color',colors{condi});
    title(['GrandRP Slope fit for all conditions, ' ' Channel ' int2str(channels(j)) ]);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    %legend([s1],'2sec','4sec','8ec','16sec','Inf','Location','Best')% COLORS DO NOT MATCH
    hold on
    filename= ['RP_grandslope_Chan_' int2str(channels(j)) '_cond_' int2str(condi), '.png'];
    saveas(k3,filename); %HOW TO SAVE ONLY LAST FIGURE FROM A LOOP??
    hold on
    hold off
end

% LEGEND NOT WORKING

%% END
%%%%%%%%%%%%%%%%%%%
