% Histograms

 ALL_TRIALS= struct('cond1',[],'cond2', [],'cond3',[], 'cond4',[], 'cond5',[]);
 
for subi= 1:22
  
    ALL_TRIALS.cond1= [ALL_TRIALS.cond1 
        pickupBehav(subi).good_resps_cond{1}'];
    ALL_TRIALS.cond2= [ALL_TRIALS.cond2 
        pickupBehav(subi).good_resps_cond{2}'];
    ALL_TRIALS.cond3= [ALL_TRIALS.cond3 
        pickupBehav(subi).good_resps_cond{3}'];
    ALL_TRIALS.cond4= [ALL_TRIALS.cond4 
        pickupBehav(subi).good_resps_cond{4}'];
    ALL_TRIALS.cond5= [ALL_TRIALS.cond5 
        pickupBehav(subi).good_resps_cond{5}'];
  
end

% ALL_TRIALS= zeros(:,5);
% ALL_TRIALS=  struct2cell(ALL_TRIALS)
% ALL_TRIALS(:,2)= ALL_TRIALS.cond2;
% ALL_TRIALS(:,3)= ALL_TRIALS.cond3;
% ALL_TRIALS(:,4)= ALL_TRIALS.cond4;
% ALL_TRIALS(:,5)= ALL_TRIALS.cond5;
%  
% ALL_TRIALS= zeros(880,5);
% for condi= 1:5
%     for subi= 1:22
%         ALL_TRIALS(:,condi)= [ALL_TRIALS(:,condi)
%             pickupBehav(subi).good_resps_cond{condi}'];
%     end
% end

% ALL_TRIALS= {};
%  
% for condi= 1:5
%     for subi= 1:22
%         ALL_TRIALS{condi}= [ALL_TRIALS{condi}
%             pickupBehav(subi).good_resps_cond{condi}'];
%     end
% end

% ALL_TRIALS= NaN(880,5);
% ALL_TRIALS= [ALL_TRIALS.cond1 ALL_TRIALS.cond2 ALL_TRIALS.cond3 ALL_TRIALS.cond4 ALL_TRIALS.cond5];


Colors= {[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.4940 0.1840 0.5560]};
Cond_names= [2 4 8 16 Inf];
Posiz1= 1:5;
Posiz2= 6:10;
figure('color','white');

for condi= 1:length(Posiz1)   
%     figure;
    ax(condi)=subplot(2,5,condi);
    subplot(ax(condi));
    title(['Condition ' num2str(Cond_names(condi)) ' sec']);
    axis([ax(condi).XLim 0 0.6]);
%     set(gcf,'color','w');
    for subi= 1:22
        
        h1(condi)=histogram(pickupBehav(subi).good_resps_cond{condi},'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',0.6); %'Normalization','count','BinMethod','auto'
%           h(condi).NumBins= 20;
%         h(condi).BinWidth= 0.5;
        %         NBins= morebins(h(condi));
        %         histogram(pickupBehav(subi).normRESPCond{condi});
        hold on;
    end
    hold on;
end   

hold on;

for condi= 1:length(Posiz2)
        %     figure;
        ax(Posiz2(condi))=subplot(2,5,Posiz2(condi));
        subplot(ax(Posiz2(condi)));
        title(['Condition ' num2str(Cond_names(condi)) ' sec']);
        axis([ax(Posiz2(condi)).XLim 0 0.35]);
        
        for subi= 1:22
            
            %         h(condi).NumBins = 15;
            %         h(condi).BinEdges = [0:6];
            h2(condi)=histogram(pickupBehav(subi).normRESPCond{condi},'BinWidth',0.25,'Normalization','probability','FaceColor',Colors{condi},'FaceAlpha',1); %'Normalization','count','BinMethod','auto'
            %     h1.NumBins= 20;
            %         h1.BinWidth= 0.1;
            %         NBins= morebins(h(condi));
            %         histogram(pickupBehav(subi).normRESPCond{condi});
            hold on;
        end
        hold on;
end

%%%%%
for subi= 1: nSubjs
    figure('color','white');
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

ax(6)=subplot(2,5,6);
subplot(ax(6));
title(['Condition 1 (' num2str(Cond_names(1)) ' sec)']);
xlabel('Waited time, normalized (sec)');
ylabel('Frequency of response');
ax(6).YLim= [0 100];
hold on;
h6=histogram(normalize(ALL_TRIALS.cond1),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{1},'FaceAlpha',0.6); 
hold on;
x1 = nanmedian(ALL_TRIALS.cond1(1));
y1=0; y2= 250;
line([x1 x1], [y1 y2],'Color','k','LineWidth',2); % median

ax(7)=subplot(2,5,7);
subplot(ax(7));
title(['Condition 2 (' num2str(Cond_names(2)) ' sec)']);
xlabel('Waited time, normalized (sec)');
ylabel('Frequency of response');
ax(7).YLim= [0 100];
hold on;
h7=histogram(normalize(ALL_TRIALS.cond2),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{2},'FaceAlpha',0.6); 
hold on;


ax(8)=subplot(2,5,8);
subplot(ax(8));
title(['Condition 3 (' num2str(Cond_names(3)) ' sec)']);
xlabel('Waited time, normalized (sec)');
ylabel('Frequency of response');
ax(8).YLim= [0 100];
hold on;
h8=histogram(normalize(ALL_TRIALS.cond3),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{3},'FaceAlpha',0.6); 
hold on;

ax(9)=subplot(2,5,9);
subplot(ax(9));
title(['Condition 4 (' num2str(Cond_names(4)) ' sec)']);
xlabel('Waited time, normalized (sec)');
ylabel('Frequency of response');
ax(9).YLim= [0 100];
hold on;
h9=histogram(normalize(ALL_TRIALS.cond4),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{4},'FaceAlpha',0.6); 
hold on;

ax(10)=subplot(2,5,10);
subplot(ax(10));
title(['Condition 5 (' num2str(Cond_names(5)) ' sec)']);
xlabel('Waited time, normalized (sec)');
ylabel('Frequency of response');
ax(10).YLim= [0 100];
hold on;
h10=histogram(normalize(ALL_TRIALS.cond5),'BinWidth',0.25,'Normalization',NM,'FaceColor',Colors{5},'FaceAlpha',0.6); 
hold on;



% test for normal distribution

[h,p,k,c] = kstest(x,'Tail','larger')
