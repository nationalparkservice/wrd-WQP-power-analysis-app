Setup
 
R packages can be built as a source package (all computers) or binary (Windows).  The npsMK R package will have a suffix of .zip or .tgz for a binary package.
 
Package deployment can be accomplished by downloading and installing the zipped package or potentially installing directly from GitHub if the user has the appropriate GitHub privileges.
 
If installing through the downloaded zipped package:
 
1. 	Open up the zipped R package file and save to your desktop.
2. 	Open a new RStudio session.
3. 	 In the RStudio console, install the package using the following line of code and the corresponding file path where the package is saved on your desktop.
 
install.packages("YOUR_FILE_PATH_HERE", type = "binary", repos = NULL)
 
For example,
 
install.packages("C:/Users/rlamont/OneDrive - DOI/Desktop/npsMK_0.0.8.tgz", type = "binary”)
 
Note the direction of the forward slashes inside the path name.
If you copy and paste the code into the console, delete and retype the quotes.
 
4. 	In the RStudio console, run the following line of code to install the remaining packages.

install.packages(c("shiny", "Kendall", "tidyverse", "dataRetrieval", "kableExtra"))
 
 
Shiny App
 
To open the app, run the following line of code in your R console.
 
source(file.path(find.package("npsMK"),"app.R"))

 
 Alternatively, open the zipped package folder and open the app.R script inside R studio. Press “Run App” in the top right corner of the script window. 

All actions are behind a “Perform Analysis” button. The user can select an NPS network, site, and parameter via the drop down menus to retrieve historical nutrient data from STORET. Data are queried from the Water Quality Portal and partitioned to simulate various sampling frequencies: annually, biannually, quarterly, and monthly. 

These different sampling frequencies are simulated using historical data to take a closer look at the relationship between sample size and power for nutrient data. The coefficient of variation (CV) is estimated for each of these frequency subsets, and the sample size (N) is recorded. Hypothetical power for each frequency scenario was calculated from a Mann-Kendall test, which is a statistical test designed to detect monotonic temporal trends. Each of the four sampling scenarios are examined through a Monte Carlo simulation using a synthetic time series with trend and a significance level of 0.05. The simulation for each frequency scenario incorporates the respective CV and a varying sample size that includes 1 to the respective N.

Varying trend slopes of 0.1, 1, and 10 are included in the simulation to allow a full visualization of the hypothetical power curve with different combinations of parameters. The user can use the slider to select an additional trend slope. 

After selecting the sampling frequency of interest, network, site, parameter, and additional trend slope, the user can press “Perform Analysis” to render plots of the hypothetical power curve. An asterisk on the X-axis corresponds to the sample size (N) of the selected sampling frequency. 
