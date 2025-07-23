clear all; clc; close all;
load struct_nomadbs_OPT1_multiUE_FullPartialCSIUMiURA64_32Wt_500R_2_45_28f500_mult
mu=2:45;
figure
plot(mu,sumrate_DBS(mu),'g','LineWidth',1.4)
hold on 
plot(mu,sumrate_CB(mu),'y','LineWidth',1.4)
hold on
plot(mu,sumrate_ZF(mu),'c','LineWidth',1.4)
hold on
plot(mu,sumrate_dbsnoma_2UE(mu),'b','LineWidth',1.4)
hold on
plot(mu,sumrate_dbsnoma_multiUE(mu),'r','LineWidth',1.4)
legend('DBS','CB','ZF','2UE-NOMA-DBS','multi UE-NOMA-DBS')
title('UMi, 28GHz, URA(2*32)')
clear all; clc; 
load struct_nomadbs_OPT1_multiUE_FullPartialCSIUMiURA64_16Wt_500R_2_45_28f500_mult
mu=2:45;
figure
plot(mu,sumrate_DBS(mu),'g*-')
hold on 
plot(mu,sumrate_CB(mu),'y*-')
hold on 
plot(mu,sumrate_ZF(mu),'c*-')
hold on
plot(mu,sumrate_dbsnoma_2UE(mu),'b*-')
hold on
plot(mu,sumrate_dbsnoma_multiUE(mu),'r*-')
legend('DBS','CB','ZF','2UE-NOMA-DBS','multi UE-NOMA-DBS')
title('UMi, 28GHz, URA(4*16)')
clear all; clc; 
load struct_nomadbs_OPT1_multiUE_FullPartialCSIUMiURA64_8Wt_500R_2_45_28f500_mult
mu=2:45;figure
plot(mu,sumrate_DBS(mu),'go-')
hold on 
plot(mu,sumrate_CB(mu),'yo-')
hold on
plot(mu,sumrate_ZF(mu),'co-')
hold on
plot(mu,sumrate_dbsnoma_2UE(mu),'bo-')
hold on
plot(mu,sumrate_dbsnoma_multiUE(mu),'ro-')
legend('DBS','CB','ZF','2UE-NOMA-DBS','multi UE-NOMA-DBS')
title('UMi, 28GHz, URA(8*8)')