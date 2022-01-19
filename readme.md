Readme file

If all the files are kept in the same directory and the working directory in the rmd file is set to the source directory, 
there should be no problems loading the functions.

Knitting the rmd file takes 10 till 15 minutes because there is a lot of code.

A bayesian regression analysis was performed and a report and the code can be found in this repository.

There is a file that contains the self-made functions that were used (functions.R), an html file with the findings, a markdown
file whith the findings, a data file and a file containing the reference that was used in the report (bibliography.bib).

In order to be able to run all the code, several pacakges must be installed. The code below can be copied and run in order to have the
necessary packages installed.
It could be possible that you get different numbers when you run the code. THis could be due to the fact
that wrong versions of packages are installed. The version I used are shown below.
install.packages('ggplot2') #(version 3.3.5)
install.packages('gridExtra') #(version 2.3)
install.packages('grid') #(version 4.0.2)
install.packages('knitr') #(version 1.37)
install.packages('kableExtra') #(version 1.3.4)
install.packages('bain') #(version 0.2.8)
install.packages('mvtnorm') #(version 1.1-3)
