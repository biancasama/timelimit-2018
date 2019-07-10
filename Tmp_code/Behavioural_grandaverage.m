%% Behav analysis script for Timelimit experiment 2.0
%% BT, created 5-6 June 2018; modified:20(?)June 2018, 21,24 August 2018.

%{
   This script will analyze behavioral data ACROSS PARTICIPANTS
   Currently the data are in /Volumes/USB DISK/TIMELIMIT_backup/Behavioural_MATLAB_output'

   THERE ARE LINES of CODE THAT NEED TO BE FIXED -- SEE BELOW

%}

%% START: initialize path
%% Load data

parent_folder='/Volumes/USB_DISK/SCIENCE/TIMELIMIT_backup/Behavioural_MATLAB_output'; %'/Users/bt_neurospin/Desktop/TimeLim_Behav'
subj_folders = dir(fullfile(parent_folder, 'Subject*'))
nSubjects= length(subj_folders);
resultsfolder = fullfile(parent_folder,'/Results');
cd(resultsfolder);

%%
RT_allsubj= [];
for subj=1:nSubjects;
    current_subject = fullfile(parent_folder, subj_folders(subj).name)
    load([current_subject,'/RTime.mat'])
    load([current_subject,'/DurationResp.mat'])
    RT_allsubj = [RT_allsubj
        RTime]; %concatenate all matrixes for single subjects (10x20)x10
end
% 
save RT_allsubj RT_allsubj 
%%
DurationResp_allsubj= [];
for subj=1:nSubjects;
    current_subject = fullfile(parent_folder, subj_folders(subj).name)
    load([current_subject,'/DurationResp.mat'])
    DurationResp_allsubj = [DurationResp_allsubj
        DurationResp]; %concatenate all matrixes for single subjects (10x20)x10
end

save DurationResp_allsubj DurationResp_allsubj
%%
conditions_all = unique(DurationResp_allsubj); %sorting
n_conditions = length(conditions_all);

save conditions_all conditions_all
%%
for condi = 1:n_conditions
    % ONLY for plot purposes
    ind_arr_cond1= find(DurationResp_allsubj == conditions_all(1));
    ind_arr_cond2= find(DurationResp_allsubj == conditions_all(2));
    ind_arr_cond3= find(DurationResp_allsubj == conditions_all(3));
    ind_arr_cond4= find(DurationResp_allsubj == conditions_all(4));
    ind_arr_cond5= find(DurationResp_allsubj == conditions_all(5));
end
    
Trlsall_1= RT_allsubj(ind_arr_cond1);
Trlsall_2= RT_allsubj(ind_arr_cond2);
Trlsall_3= RT_allsubj(ind_arr_cond3);
Trlsall_4= RT_allsubj(ind_arr_cond4);
Trlsall_5= RT_allsubj(ind_arr_cond5);

% Added 25th Oct: 
% Max Waiting times
max(Trlsall_5); max(Trlsall_4);max(Trlsall_3);max(Trlsall_2);max(Trlsall_1);
% Min 
min(Trlsall_5); min(Trlsall_4);min(Trlsall_3);min(Trlsall_2);min(Trlsall_1);

% Trials in which subject waited more than max timelimit or less than min
% timelimit
Trlsall_5(find(Trlsall_5 > 16))
Trlsall_5(find(Trlsall_5 < 2))

% Subject responded faster than 200ms
numel(Trlsall_1(find(Trlsall_1 < 0.2)))
numel(Trlsall_2(find(Trlsall_2 < 0.2)))
numel(Trlsall_3(find(Trlsall_3 < 0.2)))
numel(Trlsall_4(find(Trlsall_4 < 0.2)))
numel(Trlsall_5(find(Trlsall_5 < 0.2)))


%% Added Bootstrapped standard errors for means & medians 
condRT_allsubj = struct('mRT',[],'mdRT', [],'semRT',[]); 
for subj=1:nSubjects; % let's take 1 subject
    current_subject = fullfile(parent_folder, subj_folders(subj).name)
    load([current_subject,'/condRT.mat'])% CAVEAT (Nahuel)   
    condRT_allsubj.mRT(subj,:) = condRT.mRT;
    condRT_allsubj.mdRT(subj,:) = condRT.mdRT;
    condRT_allsubj.semRT(subj,:) = condRT.semBootRT;
    condRT_allsubj.semdRT(subj,:) = condRT.semdBootRT;

end
save condRT_allsubj condRT_allsubj
%%
conditions_all= [1 2 3 4 5];
n_conditions= length(conditions_all);

for condi= 1:n_conditions
    mean_all(condi)= mean(condRT_allsubj.mRT(:,condi))
    median_all(condi)= mean(condRT_allsubj.mdRT(:,condi))
    sem_all(condi)= mean(condRT_allsubj.semRT(:,condi))
    semd_all(condi)= mean(condRT_allsubj.semdRT(:,condi))
    mean_all= normalize(mean_all,'center','mean')
    median_all= normalize(median_all,'center','mean')
    sem_all= normalize(sem_all,'center','mean')
    semd_all= normalize(semd_all,'center','mean')
end

% save mean_all mean_all
% save median_all median_all
% save sem_all sem_all
% save semd_all semd_all
normalize(A,'norm',1)
%% 1) Plot all the histograms with sorted Response Times (in order: 2,4,8,16sec, Inf)
%FIXME
% Histogram MEDIANS
figfolder = fullfile(parent_folder,'/Figures');
cd(figfolder);

h=figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
% NOTE: can we put it in a loop?
h1= histogram(Trlsall_1,20)
% h1= histogram(condRT_allsubj.mRT(:,1),20)
hold on
x = median_all(1);
y1=0; y2=max(h1.Values);
line([x x], [y1 y2],'Color','r','LineWidth',2); %RED: mean
title('2s')
% xlim([-60 60])%??
xlabel('RT 2 sec')

subplot(2,3,2)
h2= histogram(Trlsall_2,20)
hold on
x = median_all(2);
y1=0; y2=max(h2.Values);
line([x x], [y1 y2],'Color','r','LineWidth',2); %RED: mean
title('4s')
% xlim([-60 60])
xlabel('RT 4 sec')

subplot(2,3,3)
h3= histogram(Trlsall_3,20)
hold on
x = median_all(3);
y1=0; y2=max(h3.Values);
line([x x], [y1 y2],'Color','r','LineWidth',2); %RED: mean
title('8s')
% xlim([-20 20])
xlabel('RT 8 sec')

subplot(2,3,4)
h4= histogram(Trlsall_4,20)
hold on
x = median_all(4);
y1=0; y2=max(h4.Values);
line([x x], [y1 y2],'Color','r','LineWidth',2); %RED: mean
title('16s')
% xlim([-60 60])
xlabel('RT 16 sec')

subplot(2,3,5)
h5= histogram(Trlsall_5,20)
hold on
x = median_all(5);
y1=0; y2=max(h5.Values);
line([x x], [y1 y2],'Color','r','LineWidth',2); %RED: mean
title('Inf')
% xlim([-60 60])
xlabel('RT Inf')

% subplot(2,3,6)
% h6= histogram(condRT_allsubj.mRT(:,:),20)
% hold on
% histogram(condRT_allsubj.mdRT(:,:),20)
% title('ALL')
% % xlim([-60 60])
% xlabel('RT all conds')
hold on
filename= ['RTsdistribution_ALLsubjs_.png'];
saveas(h,filename)

%% 2) Plot mean & median Response Time for all participants
%FIXME
figfolder = fullfile(parent_folder,'/Figures');
cd(figfolder);

a=2; r=2;n=5;
s = a*r.^(0:n-1); %Function rule for Recursive sequence (24 Oct)

f=figure
errorbar(s,mean_all,sem_all,'--or','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean')  % RED= mean
hold on

errorbar(s,median_all,semd_all,':sqb','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median') % BLUE= median 
set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Conditions (sec)')
ylabel('Waiting Times (sec)')
title('Mean & Median WT for ALL subjs')
legend('show','Location','best')
xlim([0 Inf]) % added 24 Oct
hold on
% filename= ['RTs_Grandavg_ALLsubjs_afterBootstrp.png'];
% saveas(f,filename);
filename= ['WaitingTimes_grandavg.png']
saveas(f,filename);

%% 3) Plot median reaction time for individuals
% FIXME
figfolder = fullfile(parent_folder,'/Figures');
cd(figfolder);

l=figure('units','normalized','outerposition',[0 0 1 1])

for subj=1:nSubjects;
    plot(1:n_conditions,condRT_allsubj.mdRT(:,:)','DisplayName', sprintf('Subject %d', nSubjects(subj)))
    legend('-DynamicLegend');
    legend('show');
    drawnow;
end
hold on
plot(1:n_conditions,median_all','r','LineWidth',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median')
set(gca,'xtick',1:5, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Conditions (sec)')
ylabel('Response Times (sec)')
legend('show','Location','northwest')
hold on
filename= ['RTs_medians_ALLsubjs_.png'];
saveas(l,filename);


%% END of the script
disp('END of the script')
clear all; close all;