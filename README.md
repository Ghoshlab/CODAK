# CODAK

The instructions on how to apply the method for any cell type abundance data set are provided in the comments with the functions (CODAK_functions.R). 






######## Instructions on how to reproduce the results ########

### Author: Pratyaydipta Rudra
### Date: 12/6/2021


1. While running any of the following codes, the working directory should be set to the folder containing this ReadMe file. 

2. Please do not make any changes to the files or directories before running the codes.

3. Go to the "codes_for_paper" folder and run all the code files in the order they are numbered. Once all codes have been run, the results and plots will be generated. See below for a description of what each code does.

4. The "functions" folder contains the different functions used for the analysis. The files can be opened to check the usage of the functions, but the files should not be changed in any way.

5. The intermediate output of the anaysis will be saved in the folder "results" and its subfolders.

6. The final plots will be saved in the folder "plots".







Description of the individual code files:

1. KvLM_binary: simulation (i)
1a. KvLM_binary_lmerperm: adding the LRT-permutation method to simulation (i)
1b. KvLM_binary - rev: adding newer methods (CODAK-ED, CODAK-BC, dcor, PERMANOVA, MiRKAT)

2. KvLM_continuous: simulation S(i)
2b. KvLM_continuous - rev: adding newer methods

3. KvLM_binary_bigsig: simulation S(ii)
3a. KvLM_binary_bigsig - rev: adding newer methods

4. KvLMC_xbin_zbin: simulation (ii)
4a. KvLMC_xbin_zbin - rev: adding newer methods

5. KvLMC_xcont_zbin: simulation (Siii)
5a. KvLMC_xcont_zbin - rev: adding newer methods

6. KvLMC_dependent: simulation (iii)
6a. KvLMC_dependent - rev: adding newer methods

7. KvLMC_dependent_bigsig: simulation S(iv)
7a. KvLMC_dependent_bigsig - rev: adding newer methods

8. KvKG_binary: simulation S(v)

9. LKvGK_x_continuous: simulation S(vi)

10. Pseudocount: simulation S(viii)

11. loosim: Simulations to compare different follow up methods

12. powercurves - rev: Create all power plots

13. effectsize_plots: Create figure 4

14. Real_data_analysis: Analysis of real data

15. equivalence: Showing equivalence of different methods

16. benchmarking: Time to run










Contributors
============

-   [Pratyaydipta Rudra](https://github.com/pratyayr)
-   [Debashis Ghosh](https://github.com/ghoshd)




















