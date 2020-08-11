# [A preliminary analysis of epidemiological trends on the Isle of Wight after the launch of the Test and Trace programme](https://www.medrxiv.org/content/10.1101/2020.07.12.20151753v1.article-info)
 
The code in this repository reproduces the analysis reported in [this paper](https://www.medrxiv.org/content/10.1101/2020.07.12.20151753v1.article-info).

Using `covid-19-cases-uk.csv` and `population_by_region.csv`, the code in `TT_analysis.html` estimates incidence of new infections and the instantaneous reproduction number R in April-June across Upper Tier Local Authorities of England.

`IoW_synthetic_control.do` uses these estimates of R (saved in `time_series_R.csv`) to perform the synthetic control analysis.