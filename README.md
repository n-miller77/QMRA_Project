# QMRA_Project
Natalie Miller
December 2025
--------------
This is a QMRA model quantifying the risks of waterborne exposure to gastrointestinal illness in the Chattahoochee River - comparing before and after precipitation events 
Exposure scenario: individuals recreating in the Chattahoochee River for 1 hour the day before a precipitation event and the day after a precipitation event. 
Pathogen of Interest: Norovirus (as I was able to find a ratio for E. coli (indicator, measured by USGS) to Norovirus, which has a dose-reponse model completed as a case study on QMRAwiki.
link: https://qmrawiki.org/case-studies/norovirus-drinking-water
--------------
Data acquisition: Data was taken over the course of a year - November 2024 to November 2025 - from a USGS monitoring station, USGS-02336000 (right at the top of Atlanta on the Chattahoochee River)
Estimated levels of E. coli concentrations (FIB), uses turbidity and other values to estimate E. coli concentrations
Estimated by regression equation, water, colonies per 100 milliliters [-1 Std Deviation Interval]
link: https://waterdata.usgs.gov/monitoring-location/USGS-02336000/#period=P365D&dataTypeId=continuous-99407-1906668952&showMedian=true&showFieldMeasurements=true
--------------
Quantitative Microbial Risk Assessment as support for bathing waters profiling (Federigi, et al) gives E. coli to pathogen ratios with a uniform distribution and a min max, from which I calculated a mean value to use as the ratio value to convert E coli to Norovirus. 
link: https://www.sciencedirect.com/science/article/abs/pii/S0025326X20304367?via%3Dihub
--------------
Note: I use a few terms interchangably in this code, I try to make it clear but just for full documentation, here are what they are: pre = dry = day before precipitation event .... after = post = wet = the day after the precipitation event 
