# HF-Trends-NHANES
This repository contains the code needed to reproduce an analysis of long-term HF trends in the US, from 1988 to 2023, using a publicly available dataset (National Health and Nutrition Examination Survey; NHANES).

NHANES data can be accessed at the following link: https://wwwn.cdc.gov/nchs/nhanes/default.aspx
Please note that you will need to download data from 1988 to 2023 (including "Continuous NHANES" and "NHANES III").

The analysis uses R, a statitistical programming language that is freely available (https://www.r-project.org/). The code herein can be used to independently reproduce our findings. Please note that you will need to modify the code in order to redirect R to the folders within which you placed NHANES files. 

- The "Analysis" file contains the code needed to reproduce the results and figures
- The "Results" file contains the code needed to reproduce (most) of the results section, with the exception of some additions after revisions.
- The "multinomial_weighted_mod.R" file contains a modified function that incorporates sampling weights.

Should you have any questions related to the file or if you run into any errors while using it, please feel free to contact me at the following email address: asu.ahmed.sayed@gmail.com
