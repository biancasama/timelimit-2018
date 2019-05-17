%% Behav analysis script for Timelimit experiment 2.0
% BT, January 2018; modified 20(?)June 2018. 3 August 2018,10,24 August
% 2018, 9th November, 10th November.

%{
   This script will analyze behavioral data for block individually saved
   Currently the data are in /Users/bt_neurospin/Desktop/TimeLim_Behav
   We will start with analyzing 1 subject   % 902 for pilot; 01 for real participant
%}

%% START: initialize path


%% Choose the subject
prompt={'Which participant do you want to look at?'};
name='Subject number';
numlines=1;
answer=inputdlg(prompt,name,numlines);
subjnum= str2double(answer);

%% Move to the right folder/path

% parent_folder= '/Volumes/LaCie/128_usb_BACKUP/Project_Timelimit/Behav_DATA'
parent_folder='/Volumes/USB_DISK/SCIENCE/TIMELIMIT_backup/Behavioural_MATLAB_output'; 
subj_folders = dir(fullfile(parent_folder, 'Subject*')); % 'Subject'%TimeLimit_S*
current_subj_folder = fullfile(parent_folder, subj_folders(subjnum).name);
cd(current_subj_folder);
datafiles = dir(fullfile(current_subj_folder, 'subj*_block*'));
n_datafiles = length(datafiles);

%% let's work without session

RTime = zeros(10,20);
DurationResp= zeros(10,20); 

for subi= i:n
for block_no= 1:n_datafiles
    
    load(fullfile(current_subj_folder, datafiles(block_no).name))
    
    n_trials= length(trials);
    for trials_no=1:n_trials
    % modified August 2018
%         if strcmp(trials(trials_no).resp_code,'onTime') & (trials(trials_no).W_time) > 0 & ~isnan(trials(trials_no).W_time); % for Timelimit 2017, added on the 11th of November
                   if strcmp(trials(trials_no).resp_code,'onTime') ; % for Timelimit 2018

            RTime(block_no, trials_no)= trials(trials_no).ReactionTime;
            DurationResp(block_no, trials_no) =trials(trials_no).DurationResp; 
        else
            RTime(block_no, trials_no)= NaN;
            DurationResp(block_no, trials_no) = NaN;
        
        end
 
    end
    
end

% Find out how many trials you discarded

BADTrials= find(isnan(RTime(:,:)))'; 
disp(['For subject ' num2str(subjnum) ' there are ' num2str(numel(BADTrials)) ' trials excluded out of ' num2str(numel(RTime(:,:))) '.']);


% Sanity check
BADCorrespCond= find(isnan(DurationResp(:,:)))';
isequal(BADTrials,BADCorrespCond)

%% ELIMINATE BAD TRIALS FROM BEHAV DATA

RTime(BADTrials)=[];
DurationResp(BADTrials)=[];

%%

% DurationResp = DurationResp(:,1) %removing the redundant part % this
% doesn't make the loop work. removed.
conditions = unique(DurationResp); %sorting
n_conditions = length(conditions);

condRT = [];
for condi = 1:5
    
    condRT.condname(condi) = conditions(condi); 
%     reshapedDurationResp = reshape(DurationResp,[1 size(DurationResp,1)*size(DurationResp,2)]);
%     reshapedRTime = reshape(RTime,[1 size(RTime,1)*size(RTime,2)]);
  
    index_array = find(DurationResp == conditions(condi)); % instead of storing the variable, save it temporally to overwrite each time to get rid of the problem of mismatch of number of elements in the structure (Zafer)
    % ONLY for plot purposes
    ind_arr_cond1= find(DurationResp == conditions(1));
    ind_arr_cond2= find(DurationResp == conditions(2));
    ind_arr_cond3= find(DurationResp == conditions(3));
    ind_arr_cond4= find(DurationResp == conditions(4));
    ind_arr_cond5= find(DurationResp == conditions(5));
    
    % GO ON with the loop [added IQR, Confidence Interval via Bootstrapping - 23rd August 2018]
    condRT.mRT(condi,:) = nanmean(RTime(index_array));% error was here: RTime(condRT.index(condi),:))extra parenthesis was treating it as a row *column instead of an index.
    condRT.mdRT(condi,:) = nanmedian(RTime(index_array));
    %added 10th Nov
    condRT.std(condi,:) = std(RTime(index_array)); 
    condRT.semRT(condi,:)= std((RTime(index_array)),'omitnan')/sqrt(length(RTime(index_array)));
%     condRT.iqr(condi,:)= iqr(RTime(index_array)); % NOTE: is it serving any purpose?
    % BOOTSTRAPPING
%     condRT.ci(condi,:) = bootci(1000,@median,RTime(index_array));
%     condRT.ci(condi,:) = bootci(1000,{@median,RTime(index_array)},'type','student','nbootstd',100);
    condRT.bootm(condi,:)= bootstrp(10000,@mean,RTime(index_array));
    condRT.semBootRT(condi,:)= std(condRT.bootm(condi,:));
    condRT.bootmd(condi,:)= bootstrp(10000,@median,RTime(index_array));
    condRT.semdBootRT(condi,:)= std(condRT.bootmd(condi,:));
    % added 25 Oct
    condRT.maxWT(condi,:)= max(RTime(index_array));
    condRT.minWT(condi,:)= min(RTime(index_array));

end

%% Save stuff here

%create folder if it doesn't already exist
[status, msg] = mkdir('Results2018');
cd(fullfile(current_subj_folder,'/Results2018'));

save BADTrials BADTrials;
save RTime RTime; save DurationResp DurationResp
save condRT condRT;

%% Plot all the histograms with sorted Response Times (in order: 2,4,8,16sec, Inf)
% Include means and medians 

% for blocki=1:n_conditions
%     
%     RTsorted(blocki,:) = RTime(condRT.index(blocki,:));
%     figure
%     histval=histogram(RTsorted(blocki,:),10)%check number of bins
%     hold on
%     x = nanmean(RTsorted(blocki,:));
%     y1=0; y2=max(histval.Values);
%     line([x x], [y1 y2],'Color','r','LineWidth',2); %RED: mean
%     hold on
%     z= nanmedian(RTsorted(blocki,:));
%     line([z z], [y1 y2],'Color','g','LineWidth',2); %GREEN: median
%     
% end

% h1= histogram(normalize(RTime(ind_arr_cond1)),'center','mean',20)

figure('units','normalized','outerposition',[0 0 1 1])
ax(1)= subplot(2,3,1)
h1= histogram(RTime(ind_arr_cond1),20)
hold on
x1 = condRT.mRT(1);
y1=0; y2=max(h1.Values);
line([x1 x1], [y1 y2],'Color','r','LineWidth',2); %RED: mean
hold on
x2 = condRT.mdRT(1);
y3=0; y4=max(h1.Values);
line([x2 x2], [y3 y4],'Color','g','LineWidth',2); %GREEN: median
title('2s')
% xlim([-5 max(RTime])
xlabel('RT 2 sec')

ax(2)= subplot(2,3,2)
h2= histogram(RTime(ind_arr_cond2),20)
hold on
x1 = condRT.mRT(2);
y1=0; y2=max(h2.Values);
line([x1 x1], [y1 y2],'Color','r','LineWidth',2); %RED: mean
hold on
x2 = condRT.mdRT(2);
y3=0; y4=max(h2.Values);
line([x2 x2], [y3 y4],'Color','g','LineWidth',2); %GREEN: median
title('4s')
% xlim([-5 60])
xlabel('RT 4 sec')

ax(3)= subplot(2,3,3);
h3= histogram(RTime(ind_arr_cond3),20);
hold on
x1 = condRT.mRT(3);
y1=0; y2=max(h3.Values);
line([x1 x1], [y1 y2],'Color','r','LineWidth',2); %RED: mean
hold on
x2 = condRT.mdRT(3);
y3=0; y4=max(h3.Values);
line([x2 x2], [y3 y4],'Color','g','LineWidth',2); %GREEN: median
title('8s')
% xlim([-20 20])
xlabel('RT 8 sec')

ax(4)= subplot(2,3,4);
h4= histogram(RTime(ind_arr_cond4),20);
hold on
x1 = condRT.mRT(4);
y1=0; y2=max(h4.Values);
line([x1 x1], [y1 y2],'Color','r','LineWidth',2); %RED: mean
hold on
x2 = condRT.mdRT(4);
y3=0; y4=max(h4.Values);
line([x2 x2], [y3 y4],'Color','g','LineWidth',2); %GREEN: median
title('16s')
% xlim([-60 60])
xlabel('RT 16 sec')

ax(5)= subplot(2,3,5);
h5= histogram(RTime(ind_arr_cond5),20);
hold on
x1 = condRT.mRT(5);
y1=0; y2=max(h5.Values);
line([x1 x1], [y1 y2],'Color','r','LineWidth',2); %RED: mean
hold on
x2 = condRT.mdRT(5);
y3=0; y4=max(h5.Values);
line([x2 x2], [y3 y4],'Color','g','LineWidth',2); %GREEN: median
title('Inf')
% xlim([-60 60])
xlabel('RT Inf')

% linkaxes(ax) % equalize axes

subplot(2,3,6)
h6= histogram(RTime,20)
hold on
x1 = condRT.mdRT(1);
y1=0; y2=max(h6.Values);
l1=line([x1 x1], [y1 y2],'Color','b','LineWidth',2); %GREEN: median
hold on
x2 = condRT.mdRT(2);
y3=0; y4=max(h6.Values);
l2=line([x2 x2], [y3 y4],'Color','r','LineWidth',2); %GREEN: median
hold on
x3 = condRT.mdRT(3);
y5=0; y6=max(h6.Values);
l3=line([x3 x3], [y5 y6],'Color','g','LineWidth',2); %GREEN: median
hold on
x4 = condRT.mdRT(4);
y7=0; y8=max(h6.Values);
l4=line([x4 x4], [y7 y8],'Color','k','LineWidth',2); %GREEN: median
hold on
x5 = condRT.mdRT(5);
y9=0; y10=max(h6.Values);
l5=line([x5 x5], [y9 y10],'Color','y','LineWidth',2); %GREEN: median
title('ALL trials')
% xlim([-60 60])
xlabel('RT all conds')
% legend(l,'2sec','3sec','5sec', '8sec','Inf') % for experiment 2017


%% Plot mean & median Response Time for individuals

% figfolder = fullfile(parent_folder,'/Figures');
% cd(figfolder);

a=2; r=2;n=5;
s = a*r.^(0:n-1); %Function rule for Recursive sequence (24 Oct)

h1=figure(1)
errorbar(s,condRT.mRT,condRT.semBootRT,'--or','LineWidth',2,'MarkerEdgeColor','r',...
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean')  % RED= mean
hold on
%Added Confidence Intervals through BOOTSTRAPPING procedure [23 August
%2018]
errorbar(s,condRT.mdRT,condRT.semdBootRT,':sqb','LineWidth',2,'MarkerEdgeColor','b',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median') % BLUE= median 
set(gca,'xtick',s, 'xticklabel',{'2s','3s','5s','8s','Inf'}) 
% set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}) 
xlabel('Time Limits Conditions (sec)')
ylabel('Waiting Times (sec)')
title(['Subj ' int2str(subjnum)])
% title(['Mean & Median RT for subj ' int2str(subjnum) ' with bootstrapped SE and Fitting']) %change it in a loop??
legend('show','Location','best')
hold on

%%
filename1= ['RTsLog_Boostr_subj_' int2str(subjnum) '.png'];
saveas(h1,filename1);

h2=figure(2);
filename2= ['Residual_fit_behav_subj' int2str(subjnum) '.png'];
saveas(h2,filename2)

%% Plot histogram from boostats

g=figure
for condi = 1:n_conditions
    histogram(condRT.boot(condi,:))
    title(['Bootstrap medians for cond ' int2str(condi)])
    hold on
    filename= ['Histogram_boostatmedians_cond ' int2str(condi) '.png'];
%     saveas(g,filename)
    hold off
end


%% End (for now)
disp(['END of the script for subj ' int2str(subjnum)])
clear all; close all;
