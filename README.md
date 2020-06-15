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

1. Run the script **[Timelimit02_ERP_average.m](https://github.com/biancasama/timelimit/blob/master/Timelimit02_ERP_average.m)** 
2. Run the script **[Timelimit02_ERP_cluster.m](https://github.com/biancasama/timelimit/blob/master/Timelimit02_ERP_cluster.m)**
3. Run the script **[Timelimit02_ERP_lateralized.m](https://github.com/biancasama/timelimit/blob/master/Timelimit02_ERP_lateralized.m)**
4. Timelimit02_ERP_plots.m
5. Timelimit02_ERP_stats.m
6. Timelimit02_ERP_statscluster.m
7. Timelimit02_ERP_statsclusterplot.m
8. Timelimit02_ERP_variability.m





---

## TRF analysis

