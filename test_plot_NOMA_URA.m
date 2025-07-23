clear all;clc;close all;
load struct_nomadbs_OPT1_multiUE_FCSIADI_FPA_FTPARMaURA64_32Wt_2000R_2_45_28f220_mult.mat
mu=3:1:45;mu_vec=3:4:45;
plot(mu,sumrate_dbsnoma_2UE_fullCSI(mu),'--c','LineWidth',2)
hold on 
plot(mu,sumrate_dbsnoma_2UE_partialCSI(mu),'k','LineWidth',2)
hold on 
plot(mu,sumrate_dbsnoma_multiUE_fullCSI(mu),'-.r','LineWidth',2)
hold on  
plot(mu,sumrate_dbsnoma_multiUE_partialCSI(mu),'g','LineWidth',2)
hold on  
plot(mu,sumrate_DBS(mu),'b','LineWidth',2)
hold on 
% load struct_nomadbs_UMiURA192_32Wt_4000R_2_45_28f500_mono
% plot(mu,sumrate_DBS(mu),'y','LineWidth',2)
% hold on 
% load struct_nomadbs_UMiURA128_32Wt_4000R_2_45_28f500_mono
% plot(mu,sumrate_DBS(mu),'m','LineWidth',2)
% hold on
% plot(mu_vec,sumrate_DBS(mu_vec),'m^','LineWidth',2)
grid on
xlabel('Number of users')
ylabel('Average system spectral efficiency [bps/Hz]')
