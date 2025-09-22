# Code and Data for OCRL

Overview

1. Directory `Simulation` contains the main functions required to perform simulations for the proposed online change rate learnig (OCRL), as well as several classical offline and online comparative methods. Specifically,

1.1 The sub-directory `OCRL` contains the functions for the OCRL algorithm.

       The R script `basic functions.R` provides the essential functions for both OCRL and the classical local polynomial regression.

       The R scripts `Mean1-st_dense_OCRL.R` and `Mean1-st_sparse_OCRL.R` estimate the first-order derivative of the mean function under dense and sparse settings.

       The R scripts `Cov-1st_C_dense_OCRL.R`, `Cov-1st_C_sparse_OCRL.R` and `Cov1-st_dense_sparse_OCRL.R` estimate the first-order derivative of the covariance function under the dense and sparse settings.

1.2 The sub-directory `Kernel` contains the functions for the classical local ploynomial regression. 

      The R scripts `Mean1-st_dense_kernel.R`, `Mean1-st_sparse_kernel.R`, `Cov-1st_C_dense_kernel`, `Cov-1st_C_sparse_kernel.R` and `Cov1-st_dense_sparse_kernel.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

1.3 The sub-directory `SGD` contains the functions for the stochastic gradient descent method.

      The R scripts  `SGDo_Mean1-st.R` and `SGDo_Cov1-st.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

1.4 The sub-directory `Spline` contains the functions for conducting simulations using the spline-based method.

      The R scripts `Spline_Mean1-st_fun.R`  and `Spline_Cov1-st_fun.R` provide the essential functions for the spline-based method.

      The R scripts `Spline_Mean1-st.R` and `Spline_Cov1-st_fun.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

1.5 The sub-directory `Wavelet` contains the functions for conducting simulations using the wavelet-based method.
      
      The R script `Waveopt_Cov1-st.R` provides the essential functions for the spline-based method.

      The R scripts `Waveopt_Mean1-st.R` and `Waveopt_Cov1-st_fun.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

2. Directory `Realdata` contains the prerocess for the two real data examples which are available in the following links:

- ***Power_Consumption*** : www.pjm.com.
- ***US_Traffic_Accident*** :  https://smoosavi.org/datasets/us_accidents.

      The R scripts `dense_mean1-st_online and batch.R`and `sparse_mean1-st_online and batch.R` provide the codes for analyzing the hourly power consumption data in megawatts.
      
      The R scripts `count_USAcc_dense_mean1-st_online and batch.R` and `count_USAcc_sparse_mean1-st_online and batch.R` provide the codes for analyzing the traffic accident data from the United States.





