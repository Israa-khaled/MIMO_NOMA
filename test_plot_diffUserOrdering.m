clear all; clc; close all;
load struct_ADNOMA_2UE_diffUserOrdering_FPARMaULA32_1Wt_600R_2_45_28f500_mult.mat
mu=2:2:45;
figure
plot(mu,SR_2UE_FPA_FCSI(mu),'g','LineWidth',1.4)
hold on 
plot(mu,SR_2UE_FPA_ADI(mu),'--k','LineWidth',1.4)
hold on 
plot(mu,SR_2UE_FPA_HW(mu),':b','LineWidth',2.2)
hold on 
plot(mu,SR_2UE_FPA_dist(mu),'-.r','LineWidth',1.4)
grid minor
xlabel('Number of users, K')
ylabel('spectral effciency [bps/Hz]')
legend('FCSI-based \zeta','ADI-based \zeta','received Power','distance')