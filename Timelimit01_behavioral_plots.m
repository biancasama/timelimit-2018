%% Plot behavioural data

%% plots

figure; bar(GAVGbehav.maxWT,'r','EdgeColor','k','LineWidth',1.5); hold on; er=errorbar(GAVGbehav.maxWT,GAVGbehav.semMaxWT); er.LineStyle = 'none'; er.LineWidth = 2.5;
title('Average maximal Waiting Times (N= 22)');
figure; bar(GAVGbehav.minWT,'EdgeColor','k','LineWidth',1.5); hold on; er=errorbar(GAVGbehav.minWT,GAVGbehav.semMinWT); er.LineStyle = 'none'; er.LineWidth = 2.5;
title('Average minimal Waiting Times (N= 22)');

% Log scale

GAVGLogbehav= [];
for condi= 1:5
    
    GAVGLogbehav.mWT(condi)= nanmean(LogBehavStats.mWT(:,condi),1);
    GAVGLogbehav.mdWT(condi)= nanmean(LogBehavStats.mdWT(:,condi),1);
    GAVGLogbehav.stdWT(condi)= nanmean(LogBehavStats.stdWT(:,condi),1);
    GAVGLogbehav.semWT(condi)= nanmean(LogBehavStats.semWT(:,condi),1);
    
    GAVGLogbehav.minWT(condi)= nanmean(LogBehavStats.minWT(:,condi),1);
    GAVGLogbehav.maxWT(condi)= nanmean(LogBehavStats.maxWT(:,condi),1);
    
end

%% Plots
% Histograms for the response times distribution across participants 
%  MEDIANS
% figfolder = fullfile(parent_folder,'/Figures');
% cd(figfolder);

h=figure('units','normalized','outerposition',[0 0 1 1])

for condi=1:5
    
    subplot(2,5,6);
    for subi= 1:length(nSubjs)
        histogram(pickupBehav(subi).RESPTIMES'); %'Normalization','count','BinMethod','auto'
        hold on;
        
    end
   
end

for condi=1:5
    
    subplot(2,5,6);
    for subi= 1:length(nSubjs)
        histogram(pickupBehav(subi).RESPTIMES'); %'Normalization','count','BinMethod','auto'
        hold on;
        
    end
   
end


x = median_all(1);
y1=0; y2=max(h1.Values);
line([x x], [y1 y2],'Color','r','LineWidth',2); %RED: mean
title('2s')
% xlim([-60 60])%??
xlabel('RT 2 sec')

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
% Plots as a function of condition
% linear (?)
a=2; r=2;n=5;
s = a*r.^(0:n-1); %Function rule for Recursive sequence (24 Oct)

% log
z= 1:5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];

hfigure=figure;
for subi= 1: length(RPsubjs)
    
    plot(s,behavStats.mdWT(subi,:),':o','Linewidth',2,...
    'MarkerSize',5);
    hold on;
    
end

% legend(h1,{'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22'},'Location','northwest');

hold on;
% function H=shadedErrorBar_seb(x,y,errBar,lineProps,transparent)
% example: shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'}); 
h2=shadedErrorBar_seb(s,behavStats.mWT(RPsubjs,:),{@mean,@sem},{'-or','LineWidth',3,'MarkerEdgeColor','r',...
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean'},1);

hold on;
h3=shadedErrorBar_seb(s,behavStats.mdWT(RPsubjs,:),{@median,@sem},{'-ob','LineWidth',3,'MarkerEdgeColor','b',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'},1);

        % errorbar(s,GAVGbehav.mWT,GAVGbehav.semWT,'-or','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean');  % RED= mean
        % hold on
        % errorbar(s,GAVGbehav.mdWT,GAVGbehav.semWT,'-sqb','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'); % BLUE= median 

% legend('SEM','mean','SEM','median','Location','northwest');
% legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','SEM','mean','SEM', 'median'},'Location','northwest');
legend({'subj 03','subj 06','subj 07','subj 08','subj 10','subj 13','subj 15','subj 17','subj 18','subj 19','subj 20','subj 21','SEM','mean','SEM', 'median'},'Location','northwest');

set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Average Waiting Times (sec)');
% xlim([0 34]); % added 24 Oct
title(['Mean & Median Waiting Times (N=' num2str(length(RPsubjs)) ')']); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for subi= 1: length(RPsubjs)
    
    
    plot(s,LogBehavStats.mdWT(subi,:),':o');
    hold on;
    
end
hold on;
%  legend(h1,{'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22'},'Location','northwest');

hold on;
% function H=shadedErrorBar_seb(x,y,errBar,lineProps,transparent)
% example: shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'}); 
h2=shadedErrorBar_seb(s,LogBehavStats.mWT(RPsubjs,:),{@mean,@sem},{'-or','LineWidth',3,'MarkerEdgeColor','r',...
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean'},1);

hold on;
h3=shadedErrorBar_seb(s,LogBehavStats.mdWT(RPsubjs,:),{@median,@sem},{'-ob','LineWidth',3,'MarkerEdgeColor','b',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'},1);

        % errorbar(s,GAVGbehav.mWT,GAVGbehav.semWT,'-or','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean');  % RED= mean
        % hold on
        % errorbar(s,GAVGbehav.mdWT,GAVGbehav.semWT,'-sqb','LineWidth',2,'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median'); % BLUE= median 

% legend('SEM','mean','SEM','median','Location','northwest');
% legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','SEM','mean','SEM', 'median'},'Location','northwest');
legend({'subj 03','subj 06','subj 07','subj 08','subj 10','subj 13','subj 15','subj 17','subj 18','subj 19','subj 20','subj 21','SEM','mean','SEM', 'median'},'Location','northwest');

set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Average Waiting Times (sec)');
% xlim([0 34]); % added 24 Oct
title(['Mean & Median Waiting Times (N=' num2str(length(RPsubjs)) ')']);

%%
figure;
for subi= 1: nSubjs
    
    plot(s,LogBehavStats.mdWT(subi,:),'o');
    hold on;
    
end
hold on;
% plot(s,GbehavStats.mWT,':or','LineWidth',2);
errorbar(s,.mWT,GAVGbehav.semWT,'-or','LineWidth',2,'MarkerEdgeColor','k',...p
    'MarkerFaceColor','red','MarkerSize',5,'DisplayName','Mean')  % RED= mean
hold on
errorbar(s,GAVGbehav.mdWT,GAVGbehav.semWT,'-sqb','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Median') % BLUE= median 

set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Waiting Times (sec)');
legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','mean','median'},'Location','northwest');
xlim([0 34]); % added 24 Oct
title('Mean & Median Waiting Times (N=22)'); 


boxplot(behavStats.mWT);

figure;
boxplot(behavStats.mdWT);
z= 1:5;
set(gca,'xtick',z, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
xlabel('Conditions (sec)');
ylabel('Waiting Times (sec)');
title('Waiting Times per condition across all participants (N=22)'); 

figure;
for condi= 1: 5
    
    plot(behavStats.stdWT(:,condi),'-o');
    hold on;
    
end

% set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf (32s)'}) ;
% % set(gca,'ytick',behavStats.stdWT(:,1), 'yticklabel',{'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22'});
% xlabel('Conditions (sec)');
% ylabel('Waiting Times (sec)');
% legend({'subj 01','subj 02', 'subj 03','subj 04','subj 05','subj 06','subj 07','subj 08','subj 09','subj 10','subj 11','subj 12','subj 13','subj 14','subj 15','subj 16','subj 17','subj 18','subj 19','subj 20','subj 21','subj 22','mean','median'},'Location','northwest');
% xlim([0 34]); % added 24 Oct
% title('Std Waiting Times (N=22)');
