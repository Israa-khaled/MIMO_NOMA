clear all; clc; close all;
load struct_nomadbs_OPT1_multiUE_FullPartialCSI_SNRRMaULA32_1Wt_50R_0_30_28f500_mult
snr=0:30;
figure
plot(snr,sumrate_dbsnoma_2UE_fullCSI,'b')
hold on 
plot(snr,sumrate_dbsnoma_2UE_partialCSI,'g')