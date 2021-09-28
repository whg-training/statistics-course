GMS Introduction to Statistics
==============================

Prerequisites
-------------
Download data and try to install the suggested software:

Suggested Software
----------------------
1. [Anaconda3 (for python3 notebooks)](https://www.anaconda.com/download/) with the following packages:
	* adjusttext - can install using `pip install adjusttext`
2. [R](https://www.r-project.org/) with the following packages:
	* FactoMineR - `install.packages("FactoMineR")` and factoextra - `install.packages("factoextra")`
	* pwr - `install.packages("pwr")`
	* XGR -
```
if(!("BiocManager" %in% rownames(installed.packages()))) install.packages("BiocManager")
BiocManager::install("XGR", dependencies=T)
```

Data
----
* Genetic distances for block 3 (contact jpwhalley if you want to download it)
* Microarray data for block 3 (contact jpwhalley if you want to download it)

Blocks
--------
1. Tests of significance and correlations
2. Dealing with uncertainty from raw experimental data and linear regression
3. Big data sets and dimension reduction methods, multiple testing and power calculations

Suggested timetable
---------
|               | Day 1        | Day 2        |
|---------------|--------------|--------------|
| 09:30 - 12:00 | Block 1      | Block 3      |
| 12:00 - 13:30 | Lunch Break  |              |
| 13:30 - 16:00 | Block 2      |     	      |

Learning Objectives
-------------------
* Those who have less experience in Statistics, confident in developing a hypothesis, choosing a test and implementing it to their data to find if significant or not.
* Those who have less experience in working with Biological datasets, get experience working with genotype data, gene expression data and raw experimental data from qPCR experiments.
* Use of confidence intervals, what to use and how to assess them in reported research.
* How to use dimension reduction methods to get an initial first look at the data
* The perils of multiple testing and the methods to test significance when doing lots of tests.

Online Resources
----------------
* [Python Tutorial](https://www.codecademy.com/learn/learn-python)
* [10 Minutes to Pandas](https://pandas.pydata.org/pandas-docs/stable/10min.html)
* [Nature Methods Points of Significance Column](https://www.nature.com/collections/qghhqm/pointsofsignificance)
* [Principal Component Methods in R: Practical Guide](http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#biplot)
* [Power analysis in R](https://www.statmethods.net/stats/power.html)
