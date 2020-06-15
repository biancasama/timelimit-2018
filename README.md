**Please, take 15 minutes to read this documentation**



---

## Behavioural analysis

You’ll start by editing this README file to learn how to edit a file in Bitbucket.

1. Run the script [**Timelimit01_behavioral.m**](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral.m). You need to have [pre-processed files](http://www.fieldtriptoolbox.org/tutorial/preprocessing_erp/) in the form '*TimeLimit_v2_Resp_subj%02d_EEG_clean_concat_rej_interp' or '*TimeLimit_2_subj%02d_EEG_clean_concat_rej_interp'.
    1a. Run my internal function **BTmy_cleandatamore** that sets up indexes for getting only the good trials and sort them based on the condition.
    1b. Load and group behavioral data, compute logaritmic transformations of the behavioural and normalize response times.
    1c. Compute *descriptive statistics*(mean, median, std, sem...).
    1d. Key outputs: *pickupBehav*, *DescriptiveStats*, containing variables: *behavStats  LogBehavStats GAVGbehav GAVGLogbehav IQR*.
    
2. If you need to control for a specific condition, run the mirror script **[Timelimit01_behavioral_control.m](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral_control.m)**.

3. Run the script **[Timelimit01_behavioral_stats.m](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral_stats.m)** to compute *inferential statistics* (Kruskalwallis test, correlations, regressions).

4. Run the script **[Timelimit01_behavioral_plots.m](https://github.com/biancasama/timelimit/blob/master/Timelimit01_behavioral_plots.m)** to make: BOX plots (documentation: [here](https://fr.mathworks.com/matlabcentral/answers/398012-adding-a-scatter-of-points-to-a-boxplot)), BAR plots, HISTOGRAMS, RAINCLOUD plots (documentation: Micah Alleh).

---

## ERP analysis





---

## TRF analysis

