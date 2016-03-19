# genetic-algorithm
## Contents
1. [Function](https://github.com/pamyo045/genetic-algorithm/blob/master/README.md#1-function)
  1. [ga_fitting.m](https://github.com/pamyo045/genetic-algorithm/blob/master/README.md#1i-ga_fittingm)
2. [Usage](https://github.com/pamyo045/genetic-algorithm/blob/master/README.md#2-usage)

***
## 1. Function
### 1.i. ga_fitting.m
This function uses MATLAB's built-in function 'ga(...)' to perform a genetic algorithm in order to find the minimum of an objective function (in this case the smallest sum of squared residual (SSR) between model parameters and experimental values within an excel file).

***
## 2. Usage
* Run ga_isotherm.m by either having ga_isotherm.m open in the MATLAB editor and pressing on the 'Run' button, or
  Enter "ga_isotherm" in MATLAB's command window and press Enter.
* Follow the instructions to select the appropriate file path, sheet, and cell range of the Excel file containing the experimental data (see 'pamyo045/genetic-algorithm/Resources/Excel Input Data Template.xlsx' for template as shown in **Fig. 1**).
  * Note: the cell range selected for the data of the Excel file must respect the template in **Fig. 1**. Make sure to respect this template style but feel free to have any number of columns (only one column of Kp values are required). For example, for **Fig. 1** the proper input for the range when prompt would be to type "A1:F22" without the quotes.
![fig1](https://github.com/pamyo045/genetic-algorithm/blob/master/Resources/Excel%20Input%20Data%20Template.png)
**Fig. 1:** Excel Input Data Template
* When prompt, decide on exporting the results to a .csv file (either already existing in the path folder of the ga_isotherm.m or enter a new name to create a new file).
  * e.g. type "results.csv" without the quotes and press Enter. This will either overwrite a file named result.csv if it already exists or create a new one if not.
