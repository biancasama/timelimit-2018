
SR = 500;
t0 = 3;
tAx = [0:2000]./SR - 3;

d = dir('subj*');
for i=1:22
%     fname = sprintf('%s/Timeseries/SpatialFilter/SF_results.mat',d(i).name);
fname= sprintf('subj%02d_SF_results_newBl', i)
    S{i} = load(fname);
end

for i=1:22
    Ind2= mean(full(S{i}.SF_timecourses_bl(S{i}.Ind2,:)));
    Ind4 = mean(full(S{i}.SF_timecourses_bl(S{i}.Ind4,:)));
    Ind8= mean(full(S{i}.SF_timecourses_bl(S{i}.Ind8,:)));
    Ind16 = mean(full(S{i}.SF_timecourses_bl(S{i}.Ind16,:)));
    IndInf = mean(full(S{i}.SF_timecourses_bl(S{i}.IndInf,:)));
    
    Ind2_bl = baseline_correct(Ind2',SR,t0,S{i}.blRange);
    Ind4_bl = baseline_correct(Ind4',SR,t0,S{i}.blRange);
    Ind8_bl = baseline_correct(Ind8',SR,t0,S{i}.blRange);
    Ind16_bl = baseline_correct(Ind16',SR,t0,S{i}.blRange);
    IndInf_bl = baseline_correct(IndInf',SR,t0,S{i}.blRange);
    
    saveTo = sprintf('peak_BL_subj%02d.mat',i);
    save(saveTo,'Ind*');
end
    
clear

for i=1:22
    fname = sprintf('peak_BL_subj%02d.mat',i);
    S(i) = load(fname);
end

for i=1:22
    MATRIX2(i,:) = S(i).Ind2_bl;
    MATRIX4(i,:) = S(i).Ind4_bl;
    MATRIX8(i,:) = S(i).Ind8_bl;
    MATRIX16(i,:) = S(i).Ind16_bl;
    MATRIXInf(i,:) = S(i).IndInf_bl;
end

% OK subjects missing

% inelegant code
% 
RPsubjs= [3  6  7  8  10  13  15  17  18  19   20   21];

    
GAVGMX2= mean(MATRIX2(RPsubjs,:));
GAVGMX4= mean(MATRIX4(RPsubjs,:));
GAVGMX8= mean(MATRIX8(RPsubjs,:));
GAVGMX16= mean(MATRIX16(RPsubjs,:));
GAVGMXInf= mean(MATRIXInf(RPsubjs,:));


GSEMMX2= sem(MATRIX2(RPsubjs,:),1);
GSEMMX4= sem(MATRIX4(RPsubjs,:));
GSEMMX8= sem(MATRIX8(RPsubjs,:));
GSEMMX16= sem(MATRIX16(RPsubjs,:));
GSEMMXInf= sem(MATRIXInf(RPsubjs,:));


mean_X2= mean(GAVGMX2(:));
mean_X4= mean(GAVGMX4(:));
mean_X8= mean(GAVGMX8(:));
mean_X16= mean(GAVGMX16(:));
mean_XInf= mean(GAVGMXInf(:));

mean_SFall= [mean_X2, mean_X4, mean_X8, mean_X16, mean_XInf];

semX2= sem(GAVGMX2(:),1);
semX4= sem(GAVGMX4(:),1);
semX8= sem(GAVGMX8(:),1);
semX16= sem(GAVGMX16(:),1);
semXInf= sem(GAVGMXInf(:),1);

sem_SFall= [semX2, semX4, semX8, semX16, semXInf];

[B,A] = butter(4,20/(SR/2));
figure
hold on
% plot(tAx,filtfilt(B,A,mean(MATRIX2(OK,:))),'linewidth',2)
% plot(tAx,filtfilt(B,A,mean(MATRIX4(OK,:))),'linewidth',2)
% plot(tAx,filtfilt(B,A,mean(MATRIX8(OK,:))),'linewidth',2)
% plot(tAx,filtfilt(B,A,mean(MATRIX16(OK,:))),'linewidth',2)

figure; hold on
plot(tAx,filtfilt(B,A,GAVGMX2(:)),'linewidth',2);
plot(tAx,filtfilt(B,A,GAVGMX4(:)),'linewidth',2);
plot(tAx,filtfilt(B,A,GAVGMX8(:)),'linewidth',2);
plot(tAx,filtfilt(B,A,GAVGMX16(:)),'linewidth',2);
plot(tAx,filtfilt(B,A,GAVGMXInf(:)),'linewidth',2);
xlabel('Time (s)');
ylabel('mean Amplitude (\muV)');
title(['RP avg N=12 with Spatial Filter']);
legend('2sec','4sec','8sec', '16sec','Inf','Location','best');

%Function rule for Recursive sequence (24 Oct)
a=2; r=2;n=5;
s = a*r.^(0:n-1); 

%% barplot
figure
bar(s,mean_SFall);
hold on
errorbar(s,mean_SFall,sem_SFall,'r.','LineWidth',2,'MarkerEdgeColor','k',...
    'MarkerFaceColor','blue','MarkerSize',5,'DisplayName','Mean');
set(gca,'xtick',s, 'xticklabel',{'2s','4s','8s','16s','Inf'}); 
xlabel('Conditions (sec)');
ylabel('RP amplitudes (\muV');
title('Spatial Filter');
