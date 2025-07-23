clear all;clc;close all;
load struct_nomadbs_OPT1_multiUE_FCSIADI_FPA_FTPARMaULA32_1Wt_5000R_2_45_28f500_mult_03
mu=3:1:45;
figure
plot(mu,sumrate_dbsnoma_2UE_FPA3(mu),':b','LineWidth',2)
hold on
load struct_nomadbs_OPT1_multiUE_FCSIADI_FPA_FTPARMaULA32_1Wt_5000R_2_45_28f500_mult_woTau 
mu=3:1:45;mu_vec=3:3:45; 
plot(mu,sumrate_dbsnoma_2UE_FPA3(mu),'--c','LineWidth',1.2)
hold on 
plot(mu,sumrate_dbsnoma_2UE_fullCSI(mu),'k','LineWidth',1.2)
hold on 
plot(mu,sumrate_dbsnoma_2UE_FPA1(mu),'-.r','LineWidth',1.2)
hold on  
plot_MarkerSize(mu,mu_vec,sumrate_dbsnoma_2UE_FTPA,'go','g')
ylim([20 60] )
grid minor
xlabel('Number of users')
ylabel('Average system spectral efficiency [bps/Hz]')
% legend('KKT-PA','FPA with \mu=0.3','FPA with \mu=0.2','FPA with \mu=0.2','FTPA')
load struct_nomadbs_OPT1_multiUE_FCSIADI_FPA_FTPARMaULA32_1Wt_5000R_2_45_28f500_mult_03
mu=3:1:45;
figure
plot(mu,sumrate_dbsnoma_2UE_FPA3_ADI(mu),':b','LineWidth',2)
hold on
load struct_nomadbs_OPT1_multiUE_FCSIADI_FPA_FTPARMaULA32_1Wt_5000R_2_45_28f500_mult_woTau 
mu=3:1:45;mu_vec=3:3:45; 
plot(mu,sumrate_dbsnoma_2UE_FPA3_ADI(mu),'k','LineWidth',1.2)
hold on 
plot(mu,sumrate_dbsnoma_2UE_partialCSI(mu),'--c','LineWidth',1.2)
hold on 
plot(mu,sumrate_dbsnoma_2UE_FPA1_ADI(mu),'-.r','LineWidth',1.2)
hold on  
plot_MarkerSize(mu,mu_vec,sumrate_dbsnoma_2UE_FTPA_ADI,'go','g')
ylim([20 60] )
grid minor
xlabel('Number of users')
ylabel('Average system spectral efficiency [bps/Hz]')
% legend('KKT-PA','FPA with \mu=0.3','FPA with \mu=0.2','FPA with \mu=0.2','FTPA')
