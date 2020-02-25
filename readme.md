# Aim
Inference of the Potts model is always problematic due to the intractable normalizing constant in its likelihood function. This project proposed a novel method to solve the normalizing constant. 
# Instruction
## folder "firstorder" is for the first order neighbourhood structure 
Below introduces the functionality of different files.
1. generatedata.m : generate data set for simulation study
2. decomnc.m   : This is the main file. Given the simulated data, MCMC is implemented. In this file, posterior samples are drawn using MCMC. 

3. composedecom.m  : Find out pixels which belong to each split. Because for each split, the spatial correlation is different.
4. showneibou.m  : This is used in composedecom.m. It is used to find out the neighbourhood of each pixel.
5. ncintegnew.m   : This is used to calculated normalizing constant according to method proposed by Peter J Green. This method is referred as TDI in our paper.
6. potts_prop.m and prop_new_potts.m  are used to generate Potts model.
7. RCoDAlike.m.   This is used to calculate the likelihood of Potts model according to our method. The method is referred as RCoDA.
8. smalldecomcover_array.m  : calculate the coverage probability.

## folder "Secondorder" is for the first order neighbourhood structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  below is instruction for folder "Secondorder"

1. generatedata.m : generate data set for simulation study
2. decomncsecond.m   : This is the main file. Given the simulated data, MCMC is implemented. In this file, posterior samples are drawn using MCMC. 

3. composedecom.m  : Find out pixels which belong to each split. Because for each split, the spatial correlation is different. This corresponds to RCoDA_M in our paper.
4. composedecom_newversion.m  : Find out pixels which belong to each split. Because for each split, the spatial correlation is different. This corresponds to RCoDA_C in our paper.
5. showneibfirst.m  : This is used in composedecom.m. It is used to find out the first order neighbourhood of each pixel.
6. showneibsecond.m  : This is used in composedecom.m. It is used to find out the second order neighbourhood of each pixel.
7. ncintegnew.m   : This is used to calculated normalizing constant according to method proposed by Peter J Green. This method is referred as TDI in our paper.
8. potts_prop.m and prop_new_potts.m  are used to generate Potts model.
7. RCoDAlike.m.   This is used to calculate the likelihood of Potts model according to our method. The method is referred as RCoDA.
8. smalldecomsecondcover_array.m  : calculate the coverage probability.
9. conschess.m : divide the Potts model into several blocks according to coding method.

#Reference
Zhu, W., & Fan, Y. (2018). A Novel Approach for Markov Random Field With Intractable Normalizing Constant on Large Lattices. Journal of Computational and Graphical Statistics, 27(1), 59-70.
