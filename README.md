**Please, take 15 minutes to read this documentation**

Follow these steps to perform analyses on M/EEG data from the Timelimit2018 experiment.

---

## Behavioural analysis

1. Run the script [**Timelimit01_behavioral.m**](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral.m). You need to have [pre-processed files](http://www.fieldtriptoolbox.org/tutorial/preprocessing_erp/) in the form '*TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp*' or '*TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp*'.
2. Run my internal function **BTmy_cleandatamore** that sets up indexes for getting only the good trials and sort them based on the condition.
3. Load and group behavioral data, compute logaritmic transformations of the behavioural and normalize response times.
4. Compute *descriptive statistics*(mean, median, std, sem...).
5. Key outputs: *pickupBehav*, *DescriptiveStats*, containing variables: *behavStats  LogBehavStats GAVGbehav GAVGLogbehav IQR*.
    
6. If you need to control for a specific condition, run the mirror script **[Timelimit01_behavioral_control.m](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral_control.m)**.

7. Run the script **[Timelimit01_behavioral_stats.m](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral_stats.m)** to compute *inferential statistics* (Kruskalwallis test, correlations, regressions).

8. Run the script **[Timelimit01_behavioral_plots.m](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral_plots.m)** to make: BOX plots (documentation: [here](https://fr.mathworks.com/matlabcentral/answers/398012-adding-a-scatter-of-points-to-a-boxplot)), BAR plots, HISTOGRAMS, RAINCLOUD plots (documentation: [Micah Alleh](https://github.com/RainCloudPlots/RainCloudPlots)).

---

## ERP analysis
Note: ER stays for event-related, TimeS for time-series.

1. Run the script **[Timelimit02_ERP_average.m](https://github.com/biancasama/timelimit/blob/master/Timelimit02_ERP_average.m)**
    This scripts computes the **movement-timelocked ERP average within subject and between subject** (grandaverage)for all channels for all time points
    within the window -3s 0s.
    After running this script you should have the following variables:
    *avg_one*: subj%02d_TimeS_one; 'one' means that all trials are mixed by condition and it's averaging across trials as well.
    *avg_trl*: subj%02d_TimeS_bytrial; 'trl' means that we are using the Fieldtrip function cfg.keeptrials  = 'yes'; to keep each trial separate (i.e. there is no averaging across trials).
    *avg_cond*: subj%02d_TimeS_cond; 'cond' means that we are sorting the trials by conditions and averaging across trials.
    *avg_condTrl*: subj%02d_TimeS_condTrl; 'condTrl' means that we are sorting the trials by conditions but using the Fieldtrip function cfg.keeptrials  = 'yes'; to keep each trial separate .
   
    *Grand_ER*: 1x5 cell array of avg 60 channels x 2001 timepoints.
    *Grand_ER_Ind*: 1x5 cell array of 22 subjects x 60 channels x 2001 timepoints (we are using the Fieldtrip function cfg.keeptrials  = 'yes'; to keep each trial separate, in this case each subject). 'Ind' stays for individual (subject).
2. Run the script **[Timelimit02_ERP_cluster.m](https://github.com/biancasama/timelimit/blob/master/Timelimit02_ERP_cluster.m)**
3. Run the script **[Timelimit02_ERP_lateralized.m](https://github.com/biancasama/timelimit/blob/master/Timelimit02_ERP_lateralized.m)**
4. Run the script **[Timelimit02_ERP_plots.m](https://github.com/biancasama/timelimit-2018/blob/master/Timelimit02_ERP_plots.m)**
5. Run the script **[Timelimit02_ERP_stats.m](https://github.com/biancasama/timelimit-2018/blob/master/Timelimit02_ERP_stats.m)**
6. Run the script **[Timelimit02_ERP_statscluster.m](https://github.com/biancasama/timelimit-2018/blob/master/Timelimit02_ERP_statscluster.m)**
7. Run the script **[Timelimit02_ERP_statsclusterplot.m](https://github.com/biancasama/timelimit-2018/blob/master/Timelimit02_ERP_statsclusterplot.m)**
8. Run the script **[Timelimit02_ERP_variability.m](https://github.com/biancasama/timelimit-2018/blob/master/Timelimit02_ERP_variability.m)**





---

## TRF analysis

