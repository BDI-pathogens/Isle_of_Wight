# [A preliminary analysis of epidemiological trends on the Isle of Wight after the launch of the Test and Trace programme](https://www.medrxiv.org/content/10.1101/2020.07.12.20151753v1.article-info)
 
The code in this repository reproduces the analysis reported in [this paper](https://www.medrxiv.org/content/10.1101/2020.07.12.20151753v1.article-info).

Using separate and combined pillars datasets, the code in `TT_analysis_1.Rmd` and displayed in [`TT_anaylsis_1.html`](https://github.com/BDI-pathogens/Isle_of_Wight/blob/master/TT_analysis_1.html)  performs a Maximum-Likelihood estimation of R over fixed time periods before and after the Test and Trace programme.

Using `pillar1_case_data.csv` and `population_by_region.csv`, the code in `TT_analysis_2.Rmd` and displayed in [`TT_anaylsis_2.html`](https://github.com/BDI-pathogens/Isle_of_Wight/blob/master/TT_analysis_2.html) estimates incidence of new infections and the instantaneous reproduction number R in April-June across Upper Tier Local Authorities of England.

The code in `IoW_synthetic_control.do` uses these estimates of R (saved in `time_series_R.csv`) to perform the synthetic control analysis.