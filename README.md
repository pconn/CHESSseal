# CHESSseal
Analysis of ice-associated seal abundance from 2016 aerial surveys in the Chukchi Sea.  This repository is associated with the paper

"Abundance and distribution of ringed and bearded seals in the Chukchi Sea: a reference for future trends"  by P.L. Boveng, V.I. Chernook, E.E. Moreland, P.B. Conn, I.S. Trukhanova, et al.  (to be published in Arctic Science sometime in 2025)

Analysis is conducted with the R script "run_CHESS_BovengEtAl2025.R" which 
requires compiling the TMB templated-C++ file, "CHESS_BovengEtAl2025cpp.cpp"
The script also relies on loading "CHESS_data_Nov2021.RData." 

In order to run the script, the user will need to download the files into their directory of choice, and to update the directory structures in "run_CHESS_BovengEtAl2025.R" to reflect this change.

An example of how bootstrapping for variance estimation was conducted is provided in "run_CHESS_BovengEtAl2025_boot.R".  In practice, bootstrapping was conducted using four versions of this file, with different random starting seeds and 150 bootstrap iterations per instance.

For those interested in raw data, please see:

Paul Conn, Irina Trukhanova, Stacie Hardy, Peter Boveng, Vladimir Burkanov, & Eric Regehr. (2021). Detections of seals, polar bears, and polar bear tracks obtained during aerial surveys over the Chukchi Sea during April and May of 2016. Research Workspace. 10.24431/rw1k588

The connection between the raw data and the R workspace used here is as follows:

Upon loading CHESS_data_nov21.RData, seal counts for US and Russian surveys are available in the `CHESS_data$C_us` and `CHESS_data$C_rus` objects; each count is associated with a specific grid cell and day surveyed.  These are held in the `CHESS_data$S_us_i` and `CHESS_data$S_rus_i` objects, where the integer value is equal to (day_i-1)*n_s+grid_i. The number of grid cells is n_s=1354, so for example if grid cell 600 was visited on day 2 of the survey, the value of S would be 1354+600 = 1954.
