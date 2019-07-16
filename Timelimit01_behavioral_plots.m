%% Plot behavioural data
%=========================================================================%
% AUTHOR: Bianca Trovo (bianca.trovo@alumni.unitn.it)
% DATE: created on July 2019
% EXPERIMENT: Timelimit_2018

%{

    SCOPE: 
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

%% More specific paths (maybe set this in the start script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

behavioral_folder= [results_Path, '/Behaviour']; % it can be also current_subj_folder
if ~exist(fullfile(behavioral_folder)); mkdir(fullfile(behavioral_folder)); end;

    
%% Load behavioral data 

cd(behavioral_folder);
load 'DescriptiveStats'; load 'pickupBehav';

%% BOX-plots

% documentation: https://fr.mathworks.com/matlabcentral/answers/398012-adding-a-scatter-of-points-to-a-boxplot

% better keep everything in log scale
myData= LogBehavStats.mdWT; % behavStats.mdWT or LogBehavStats.mdWT
figure;
H= boxplot(myData,'Notch','on','Labels',{'2s','4s','8s','16s','Inf'},'Whisker',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filled box version
% get(H,'tag');
% set(H,'Color','k','LineWidth',2);
Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]}; % rainbow scale of colors
% h= findobj(gca,'Tag','Box');
% for i=1:length(h)
%     patch(get(h(i),'XData'),get(h(i),'YData'),'Colors',Colors(h),'FaceAlpha',0.5);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% empty box but working version
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'g');
% Change the boxplot color from blue to green
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
box1 = a(7);   % The 7th object is the first box
set(a, 'Color', 'b');   % Set the color of the first box to green
hold on;
% add somerandomness 
x1=ones(length(myData(:,1))).*(1+(rand(length(myData(:,1))-0.5)/5));
x2=ones(length(myData(:,2))).*(1+(rand(length(myData(:,2)-0.5)/10));
x3=ones(length(myData(:,3))).*(1+(rand(length(myData(:,3))-0.5)/15));
x4=ones(length(myData(:,4))).*(1+(rand(length(myData(:,4))-0.5)/20));
x5=ones(length(myData(:,5))).*(1+(rand(length(myData(:,5))-0.5)/25));

x1=ones(length(myData(:,1)));
x2=ones(length(myData(:,2)));
x3=ones(length(myData(:,3)));
x4=ones(length(myData(:,4)));
x5=ones(length(myData(:,5)));

f1=scatter(x1(:,1),myData(:,1),'k','filled');f1.MarkerFaceAlpha = 0.4;hold on 
f2=scatter(x2(:,2).*2,myData(:,2),'k','filled');f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on;
f3=scatter(x3(:,3).*3,myData(:,3),'k','filled');f3.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on;
f4=scatter(x4(:,4).*4,myData(:,4),'k','filled');f4.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on;
f5=scatter(x5(:,5).*5,myData(:,5),'k','filled');f5.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on;

ylim([0,max(myData(:,5))+std(myData(:,5))]); % not sure if it's working right now...
xlabel('Conditions (sec)','FontSize',34);
ylabel('Waiting Times (log scale)','FontSize',34);
title('Waiting Times per condition across all participants (N=22)','FontSize',34); 


%% BAR plots (not used for the paper)
% Used SEM from mean and not medians but plotting medians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% median waiting times
f=figure;
ba=bar(GAVGbehav.mdWT,'EdgeColor','k','LineWidth',1.5 );
hold on; er=errorbar(GAVGbehav.mdWT,GAVGbehav.semWT); 
er.LineStyle = 'none'; er.LineWidth = 2.5;
xlabel('Timelimit conditions','Fontsize',34);
set(gca,'xtick',1:5, 'xticklabel',{'2s','4s','8s','16s','Inf'},'Fontsize',34);
ylabel('Waiting times (sec)','Fontsize',34);
title('Median Waiting Times (N= 22)','FontSize', 34);

filename= ['Barplot_WT_Timelimit2018.png'];
saveas(f,filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% average maximal waiting times 
figure; bar(GAVGbehav.maxWT,'EdgeColor','k','LineWidth',1.5); hold on; er=errorbar(GAVGbehav.maxWT,GAVGbehav.semMaxWT); er.LineStyle = 'none'; er.LineWidth = 2.5;
title('Average maximal Waiting Times (N= 22)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% average minimal waiting times 
figure; bar(GAVGbehav.minWT,'EdgeColor','k','LineWidth',1.5); hold on; er=errorbar(GAVGbehav.minWT,GAVGbehav.semMinWT); er.LineStyle = 'none'; er.LineWidth = 2.5;
title('Average minimal Waiting Times (N= 22)');



%% HISTOGRAMS
% Plots in order fo use 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters

Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]}; % rainbow scale of colors
Cond_names= [2 4 8 16 Inf]; Cond_values= [2 4 8 16 36];
Posiz1= 1:5;
Posiz2= 6:10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ALL TRIALS pooled from all subjects + try a fitting (in comment)
figure('color','white');
for condi= 1:length(Posiz1)   
    %     figure;
    ax(condi)=subplot(2,5,condi);
    subplot(ax(condi));
    
    for subi= 1:22
        
        mydata= pickupBehav(subi).good_resps_cond{condi};
        
        h1(condi)=histogram(mydata,'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.6); %'Normalization','count','BinMethod','auto'
%                   h1(condi).NumBins= 20;
%                 h1(condi).BinWidth= 0.5;
        %         NBins= morebins(h(condi));
        %         histogram(pickupBehav(subi).normRESPCond{condi});
        xlabel('Waiting times (sec)','FontSize',18);
        ylabel('Probability','FontSize',18);
        title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',18);
        set(ax(condi),'fontsize', 18);
        ax(condi).XLim= [ax(condi).XLim(1) Cond_values(condi)];
        if condi== 1
            ax(condi).YLim== [ax(condi).YLim(1) 0.7];
        end
        if condi== 1
            xticks([0:1:2]);
        elseif condi== 2
            xticks([0:2:4]);
        elseif condi== 3
            xticks([0:4:8]);
        elseif condi== 4
            xticks([0:4:16]);
        else condi== 5
            xticks([0:8:36]);
        end
    
        hold on;
        
    end
    
%         dataforfit= pickupBehav(22).good_resps_cond{condi};
%         distribution = 'Lognormal';
%         pd = fitdist(dataforfit',distribution)
%         x_values = linspace(0, median(pickupBehav(subi).good_resps_cond{condi}), 100);
%         y = pdf(pd,x_values);
%         % subplot(2, 2, 3);
%         hold on;
%         plot(x_values,y,'LineWidth',2)
      
    hold on;
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % histogram version on group level but with no adjusted bins (see below) +
% % fitting attempt
% figure('color','white');
% for condi= 1:5
%     %figure('color','white');
%     ax(condi)=subplot(2,5,condi);
%     subplot(ax(condi));
% %     axis([ax(condi).XLim 0 3]);
%     
%     myData= behavStats.mdWT(:,condi);
%     mycolors= Colors{condi};
%     
%     h3(condi)=histogram(myData,'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.7);
%     
%     xlabel('Waiting times (sec)','FontSize',18);
%     ylabel('Probability','FontSize',18);
%     title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',18);
%     set(ax(condi),'fontsize', 18);
%     ax(condi).XLim= [ax(condi).XLim(1) Cond_values(condi)];
% %     if condi== 1
% %         ax(condi).YLim== [ax(condi).YLim(1) 0.7];
% %     end
% %     if condi== 1
% %         xticks([0:1:max(behavStats.mdWT(:,1))]);
% %     elseif condi== 2
% %         xticks([0:2:max(behavStats.mdWT(:,2))]);
% %     elseif condi== 3
% %         xticks([0:4:max(behavStats.mdWT(:,3))]);
% %     elseif condi== 4
% %         xticks([0:4:max(behavStats.mdWT(:,4))]);
% %     else condi== 5
% %         xticks([0:8:max(behavStats.mdWT(:,5))]);
% %     end
% 
%     % fitting not working
% %     distribution = 'Lognormal';
% %     pd = fitdist(myData,distribution)
% %     x_values = linspace(0, max(myData),100);
% %     y = pdf(pd,x_values);
% %     % subplot(2, 2, 3);
% %     hold on;
% %     plot(x_values,y,'LineWidth',2);
% %    
%     box off;
%     hold on;
%     
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on group level (median values) - we won't probably use it
figure('color','white');
    for condi= 1:length(Posiz1)   
%     figure('color','white');
        ax(condi)=subplot(2,5,condi);
        ax(condi).FontSize=34;
        subplot(ax(condi));
        xlabel('Waiting times (sec)','FontSize',34);
        ylabel('Probability','FontSize',34);
        title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',34);
        axis([ax(condi).XLim 0 0.6]);
 
            h1(condi)=histogram(behavStats.mdWT(:,condi),'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.6); %'Normalization','count','BinMethod','auto'
            h1(condi).NumBins= 10;
            
        box off;
        hold on;
    end 
    
hold on;
% 
% % second line of histograms of normalized values: replace Indiv subjects
% % with gavg
% for condi= 1:length(Posiz2)
%         %     figure;
%         ax(Posiz2(condi))=subplot(2,5,Posiz2(condi));
%         subplot(ax(Posiz2(condi)));
%         title(['Condition ' num2str(Cond_names(condi)) ' sec']);
%         axis([ax(Posiz2(condi)).XLim 0 0.35]);
%         
%         for subi= 1:22
%             
%             %         h(condi).NumBins = 15;
%             %         h(condi).BinEdges = [0:6];
%             h2(condi)=histogram(pickupBehav(subi).normRESPCond{condi},'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.6); %'Normalization','count','BinMethod','auto'
%             h2.NumBins= 10;
%             %         h1.BinWidth= 0.1;
%             %         NBins= morebins(h(condi));
%             %         histogram(pickupBehav(subi).normRESPCond{condi});
%             hold on;
%         end
%         hold on;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot an histogram for each subject separately (adjustfor row of subplot:
% now it's just plotting one at the time...)

for subi= 1: nSubjs
    %figure('color','white');
    NM= 'count';
    ax(1)=subplot(1,5,1);
    subplot(ax(1));
    title(['Condition 1 (' num2str(Cond_names(1)) ' sec)']);
    xlabel('Waited time (sec)');
    set(gca,'xtick',1:3, 'xticklabel',{'1s','2s'}) ;
    ylabel('Frequency of response');
    ax(1).XLim= [0 2];
    % ax(1).YLim= [0 200];
    % axis([ax(1).XLim 0 250]);
    % ax1=subplot(2,5,1);
    % % subplot(ax);
    % t1=title(['Condition ' num2str(Cond_names(1)) ' sec']);
    % axis([ax(1).XLim 0 0.6]);
    hold on;
    h1(subi)=histogram(pickupBehav(subi).good_resps_cond{1}','NumBins',40,'Normalization','probability','FaceColor',Colors{1},'FaceAlpha',0.6);
    hold on;
    % x1 = nanmedian(ALL_TRIALS.cond1);
    % x1= GbehavStats.mdWT(1);
    % y1=0; y2= 250;
    % line([x1 x1], [y1 y2],'Color','k','LineWidth',2); % median
    
    ax(2)=subplot(1,5,2);
    subplot(ax(2));
    title(['Condition 2 (' num2str(Cond_names(2)) ' sec)']);
    xlabel('Waited time (sec)');
    ylabel('Frequency of response');
    ax(2).XLim= [0 4];
    % ax(2).YLim= [0 200];
    % ax2=subplot(2,5,2);
    % % subplot(ax(1));
    % t2=title(['Condition ' num2str(Cond_names(2)) ' sec']);
    % axis([ax(1).XLim 0 0.6]);
    hold on;
    h2=histogram(pickupBehav(subi).good_resps_cond{2}','NumBins',40,'Normalization','probability','FaceColor',Colors{2},'FaceAlpha',0.6);
    hold on;
    % x2 = nanmedian(ALL_TRIALS.cond2(2));
    % x2= GbehavStats.mdWT(2);
    % y1=0; y2= 250;
    % line([x2 x2], [y1 y2],'Color','k','LineWidth',2); % median
    
    ax(3)=subplot(1,5,3);
    subplot(ax(3));
    title(['Condition 3 (' num2str(Cond_names(3)) ' sec)']);
    xlabel('Waited time (sec)');
    ylabel('Frequency of response');
    ax(3).XLim= [0 8];
    % ax(3).YLim= [0 200];
    hold on;
    h3=histogram(pickupBehav(subi).good_resps_cond{3}','NumBins',40,'Normalization','probability','FaceColor',Colors{3},'FaceAlpha',0.6);
    hold on;
    % x3 = nanmedian(ALL_TRIALS.cond3(3));
    % x3= GbehavStats.mdWT(3);
    % y1=0; y2= 250;
    % line([x3 x3], [y1 y2],'Color','k','LineWidth',2); % median
    
    ax(4)=subplot(1,5,4);
    subplot(ax(4));
    title(['Condition 4 (' num2str(Cond_names(4)) ' sec)']);
    xlabel('Waited time (sec)');
    ylabel('Frequency of response');
    ax(4).XLim= [0 16];
    % ax(4).YLim= [0 200];
    hold on;
    h4=histogram(pickupBehav(subi).good_resps_cond{4}','BinWidth',0.25,'Normalization','probability','FaceColor',Colors{4},'FaceAlpha',0.6);
    hold on;
    % x4 = nanmedian(ALL_TRIALS.cond4(4));
    % x4= GbehavStats.mdWT(4);
    % y1=0; y2= 250;
    % line([x4 x4], [y1 y2],'Color','k','LineWidth',2); % median
    
    
    ax(5)=subplot(1,5,5);
    subplot(ax(5));
    title(['Condition 5 (' num2str(Cond_names(5)) ' sec)']);
    xlabel('Waited time (sec)');
    ylabel('Frequency of response');
    ax(5).XLim= [0 max(pickupBehav(subi).good_resps_cond{5}')];
    % ax(5).YLim= [0 200];
    hold on;
    h5=histogram(pickupBehav(subi).good_resps_cond{5}','NumBins',40,'Normalization','probability','FaceColor',Colors{5},'FaceAlpha',0.6);
    hold on;
    % x5 = nanmedian(ALL_TRIALS.cond5(5));
    % x5= GbehavStats.mdWT(5);
    % y1=0; y2= 250;
    % line([x5 x5], [y1 y2],'Color','k','LineWidth',2); % median
end

% % to plot all trials pooled together from all subjects (is this correct??)
% 
% ALL_TRIALS= struct('cond1',[],'cond2', [],'cond3',[], 'cond4',[], 'cond5',[]);
%  
% for subi= 1:22
%   
%     ALL_TRIALS.cond1= [ALL_TRIALS.cond1 
%         pickupBehav(subi).good_resps_cond{1}'];
%     ALL_TRIALS.cond2= [ALL_TRIALS.cond2 
%         pickupBehav(subi).good_resps_cond{2}'];
%     ALL_TRIALS.cond3= [ALL_TRIALS.cond3 
%         pickupBehav(subi).good_resps_cond{3}'];
%     ALL_TRIALS.cond4= [ALL_TRIALS.cond4 
%         pickupBehav(subi).good_resps_cond{4}'];
%     ALL_TRIALS.cond5= [ALL_TRIALS.cond5 
%         pickupBehav(subi).good_resps_cond{5}'];
%   
% end
% 
% figure('color','white');
% % here normaliwed values
% ax(6)=subplot(2,5,6);
% subplot(ax(6));
% title(['Condition 1 (' num2str(Cond_names(1)) ' sec)']);
% xlabel('Waited time, normalized (sec)');
% ylabel('Frequency of response');
% ax(6).YLim= [0 100];
% hold on;
% h6=histogram(normalize(ALL_TRIALS.cond1),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{1},'FaceAlpha',0.6); 
% hold on;
% x1 = nanmedian(ALL_TRIALS.cond1(1));
% y1=0; y2= 250;
% line([x1 x1], [y1 y2],'Color','k','LineWidth',2); % median
% 
% ax(7)=subplot(2,5,7);
% subplot(ax(7));
% title(['Condition 2 (' num2str(Cond_names(2)) ' sec)']);
% xlabel('Waited time, normalized (sec)');
% ylabel('Frequency of response');
% ax(7).YLim= [0 100];
% hold on;
% h7=histogram(normalize(ALL_TRIALS.cond2),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{2},'FaceAlpha',0.6); 
% hold on;
% 
% 
% ax(8)=subplot(2,5,8);
% subplot(ax(8));
% title(['Condition 3 (' num2str(Cond_names(3)) ' sec)']);
% xlabel('Waited time, normalized (sec)');
% ylabel('Frequency of response');
% ax(8).YLim= [0 100];
% hold on;
% h8=histogram(normalize(ALL_TRIALS.cond3),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{3},'FaceAlpha',0.6); 
% hold on;
% 
% ax(9)=subplot(2,5,9);
% subplot(ax(9));
% title(['Condition 4 (' num2str(Cond_names(4)) ' sec)']);
% xlabel('Waited time, normalized (sec)');
% ylabel('Frequency of response');
% ax(9).YLim= [0 100];
% hold on;
% h9=histogram(normalize(ALL_TRIALS.cond4),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{4},'FaceAlpha',0.6); 
% hold on;
% 
% ax(10)=subplot(2,5,10);
% subplot(ax(10));
% title(['Condition 5 (' num2str(Cond_names(5)) ' sec)']);
% xlabel('Waited time, normalized (sec)');
% ylabel('Frequency of response');
% ax(10).YLim= [0 100];
% hold on;
% h10=histogram(normalize(ALL_TRIALS.cond5),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{5},'FaceAlpha',0.6); 
% hold on;


%% Plots as a function of condition
% linear (?)
a=2; r=2;n=5;
s = a*r.^(0:n-1); %Function rule for Recursive sequence (24 Oct)

% log
z= 1:5;
RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% non-log version
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

% semi-log version
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


%% Standard deviations
mycolormap= [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; 0 0.4470 0.7410; 0.4940 0.1840 0.5560]; 
myData= behavStats.stdWT; % behavStats.stdWT or LogBehavStats.mdWT

f1=figure('color','white');
for condi= 1: 5
    
    p=plot(myData(:,condi),'o','Color',mycolormap(condi,:),'LineWidth',1,'MarkerSize',15,'MarkerFaceColor',mycolormap(condi,:));
   
    hold on;
    
end
% Create rectangle
annotation(f1,'rectangle',...
    [0.150498997995992 0.11756061719324 0.0202424849699399 0.120499632623071],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',3);

% Create rectangle
annotation(f1,'rectangle',...
    [0.24348496993988 0.11756061719324 0.0202424849699399 0.111682586333578],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',3);

% Create rectangle
annotation(f1,'rectangle',...
    [0.367332665330663 0.11756061719324 0.0202424849699399 0.22556943423953],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',3);

% Create rectangle
annotation(f1,'rectangle',...
    [0.429456913827656 0.11756061719324 0.0202424849699399 0.22556943423953],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',3);

% Create rectangle
annotation(f1,'rectangle',...
    [0.274346693386775 0.11756061719324 0.0202424849699399 0.0830271858927261],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',3);

% Create rectangle
annotation(f1,'rectangle',...
    [0.801801603206414 0.11756061719324 0.0202424849699399 0.093313739897134],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',3);

% Create rectangle
annotation(f1,'rectangle',...
    [0.615428857715432 0.11756061719324 0.0202424849699399 0.0815576781778102],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',3);

% set(gca,'ytick',1:5, 'yticklabel',{'2s','4s','8s','16s','Inf'}) ;
box off;
set(gca,'xtick',1:22, 'xticklabel',{'01','02', '03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22'},'FontSize',28);
xlabel('Subjects (numbers)');
ylabel('Waiting Times std (sec)');
l=legend({'2s','4s','8s','16s','Inf','Non-Timers'},'FontSize',34);
% legend boxoff;
title('Std Waiting Times (N=22)','FontSize',34);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CORRELATION: Is the Standard deviation scaled with the Mean?


X= behavStats.mdWT;
Y= behavStats.stdWT;
sz = 15;
f3=figure('color','white');
for condi= 1: 5
    
%     c= mycolormap(condi,:);
    plot(X(:,condi),Y(:,condi),'o','Color',mycolormap(condi,:),'LineWidth',1,'MarkerSize',15,'MarkerFaceColor',mycolormap(condi,:))
   
    hold on;
    
end
box off;
set(gca,'FontSize',28);
xlabel('Mean Waiting Times (sec)','FontSize',34);
ylabel('Waiting Times std(sec)','FontSize',34);
title('Std is scaled with the mean','FontSize',34);
l=legend({'2s','4s','8s','16s','Inf'},'FontSize',34);
str={'R= 0.88'};
t= text(12,3,str,'FontSize',34);

%% RAINCLOUDplots

% documentation: Micah Alleh

% --------------------------
% Example: h = raincloud('X', myData, 'box_on', 1, 'color', [0.5 0.5 0.5]) 
% Parameters 
Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]};
Cond_names= [2 4 8 16 Inf];
x_limits= {[0 2],[0 4],[0 8],[0 16],[0 Inf]};
ylims= [
%'box_dodge' ,1,'box_dodge_amount',.15,'dot_dodge_amount',.35,'box_col_match',1

%% Medians

figure('color','white');
for condi= 1:5
    %figure('color','white');
    ax(condi)=subplot(2,5,condi);
    subplot(ax(condi));
%     axis([ax(condi).XLim 0 3]);
    
    myData= behavStats.mdWT(:,condi);
    mycolors= Colors{condi};
    h=raincloud_plot('X',myData,'color',mycolors,'alpha',.7,'line_width',1,'box_on', 1); 
    
    title(['Condition ' num2str(Cond_names(condi)) ' sec'],'FontSize',20);
    set(ax(condi),'fontsize', 20);
    xlabel('Waiting times (sec)', 'fontsize', 18);
%     xticks([0:.5:2]);
    box off;
    hold on;
end
    

% all trials
    for condi= 1:length(Posiz1)
        %     figure;
        ax(condi)=subplot(2,5,condi);
        subplot(ax(condi));
        title(['Condition ' num2str(Cond_names(condi)) ' sec']);
        axis([ax(condi).XLim 0 0.6]);
        %     set(gcf,'color','w');
        figure();
        for subi= 1:22
            myData= pickupBehav(subi).good_resps_cond{condi};
            mycolors= Colors{condi};
            h=raincloud_plot('X',myData,'color',mycolors,'FaceAlpha',0.6); %

            hold on;
        end
        hold on;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END
