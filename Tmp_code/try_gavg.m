SR = 500;
t0 = 3;
tAx = [0:2000]./SR - 3;

d = dir('subj*');

for i=1:22
%     fname = sprintf('%s/Timeseries/SpatialFilter/SF_results.mat',d(i).name);
fname= sprintf('subj%02d_SF_results_newBl', i)
    S{i} = load(fname);
   
% across trials 5up tp 40) and across data points (2001)?? 
mean(full(SF_timecourses_bl(Ind2,2)))
mean(full(SF_timecourses_bl(Ind4,2)))
mean(full(SF_timecourses_bl(Ind8,2)))
mean(full(SF_timecourses_bl(Ind16,2)))
mean(full(SF_timecourses_bl(IndInf,2)))

h=figure; plot(tAx,mean(full(SF_timecourses_bl(IndInf,:))),'linewidth',2)
hold on
plot(tAx,mean(full(SF_timecourses_bl(Ind16,:))),'linewidth',2)
hold on
plot(tAx,mean(full(SF_timecourses_bl(Ind8,:))),'linewidth',2)
hold on
plot(tAx,mean(full(SF_timecourses_bl(Ind4,:))),'linewidth',2)
hold on
plot(tAx,mean(full(SF_timecourses_bl(Ind2,:))),'linewidth',2)

xlabel('Time (s)');
ylabel('mean Amplitude (\muV)');
title(['Subj ' num2str(subjnum) ', RP avg with Spatial Filter']);
legend('Inf', '16sec', '8sec', '4sec', '2sec','Location','SouthWest');
% legend('Infinity condition', 'Location','NorthWest');

%  save figure for further comparisons
filename= [sprintf('Readiness_Potential_SF_subj%02d', subjnum) '.png'];
cd(figures_Path);
saveas(h,filename);
