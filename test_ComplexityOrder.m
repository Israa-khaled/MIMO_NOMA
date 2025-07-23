clear all;clc;
addpath TempFct2UE mUE
rmpath nyusimv1.6.1
scn='mult';
if scn == 'mono'
  addpath matlab_code_NYUSIM_monotraject
  rmpath matlab_code_NYUSIM_multitraject  
elseif scn == 'mult'
  rmpath matlab_code_NYUSIM_monotraject
  addpath matlab_code_NYUSIM_multitraject
end
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
load inputParameters 
f=28;sceType = 'RMa';Pe=1;Pmin=1e-3;
Nt=32;TxArrayType='ULA';Wt=1;
MU=10;% the users belongs to the same side of BS, i.e., 0<theta<180.
r=1;R=100;beta_seuil=0.5;TXPower_mw=Pe*1e3;TXPower=10*log10(TXPower_mw);NP=noise_power;
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
while (r<=R)
[CIR_MIMO_Struct,powerSpectrum,H_total1,H,AOD_LobePowerSpectrum,powerspectrum_struct,Pr_dBm,TR_dist] = getH_MIMO_n_1(f,RFBW,sceType,envType,TXPower,Nt,Nr,MU,dmin,dmax,h_BS,TxArrayType,Wt,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX);  
% DBS preocding matrix
[azimuthAOD_USER,eleAOD_USER,a_BS,W_DBSn]= getWdbs(powerspectrum_struct,Nt,Wt,TxArrayType,MU,dTxAnt);  
HW_DBS=H_total1*W_DBSn;
[SINR_user_DBS,SINR_user_CB,SINR_user_ZF,rate_user_DBS,rate_user_CB,rate_user_ZF] =getsumrateULAURA_Channelisgiven(a_BS,H_total1,MU,noise_power);
sumrate_DBS1(r)=sum(rate_user_DBS);
sumrate_CB1(r)=sum(rate_user_CB);
sumrate_ZF1(r)=sum(rate_user_ZF);
%% User Clustering Algorithm
% Select the set of groups of 2-UEs who have a spatial interference greater than \beta_0 
[Grp2b2Interf,beta_n_allUser,beta_n_seuil,beta_Grp2b2Interf]= UserbetaMatrix_beta0(a_BS,beta_seuil,MU,Nt);
%Multi-UE clustering:
[Clusters_array_multiUE,nb_2UEinCluster_multiUE,nb_clusters_multiUE,FFF_real(r)]=get_UserClustering_multiUE(MU,Grp2b2Interf,beta_Grp2b2Interf);
%2-UE clustering:
[Clusters_array_2UE,nb_2UEinCluster_2UE,nb_clusters_2UE,FFF_2UE_real(r)]=get_UserClustering_2UE(beta_n_seuil,MU);
%%
C_mltUE_1=0;CSI_mltUE_1=0;C1_mltUE_1=0;
for mm=1:nb_2UEinCluster_multiUE
    C_mltUE_1=C_mltUE_1+length(Clusters_array_multiUE{mm});
    CSI_mltUE_1=CSI_mltUE_1+length(Clusters_array_multiUE{mm})-1;
    C1_mltUE_1=C1_mltUE_1+fix(length(Clusters_array_multiUE{mm})/2);
end
C_mltUE_real(r)=C_mltUE_1;
C1_mltUE_real(r)=C1_mltUE_1;
CSI_mltUE_real(r)=CSI_mltUE_1;

nb_2UEinCluster_2UE_real(r)=nb_2UEinCluster_2UE;
nb_2UEinCluster_multiUE_real(r)=nb_2UEinCluster_multiUE;

C_2UE_real(r)=2*nb_2UEinCluster_2UE;
CSI_2UE_real(r)=nb_2UEinCluster_2UE;
FFF_2UE_realOrg(r)=(MU-1)*MU*nb_2UEinCluster_2UE;
r=r+1;
end
FFF=sum(FFF_real)/R
C1_mltUE=sum(C1_mltUE_real)/R
C_mltUE=sum(C_mltUE_real)/R
CSI_mltUE_real=sum(CSI_mltUE_real)/R

C2=sum(nb_2UEinCluster_2UE_real)/R
Cm=sum(nb_2UEinCluster_multiUE_real)/R

FFF_2UEOrg=sum(FFF_2UE_realOrg)/R
FFF_2UE=sum(FFF_2UE_real)/R
C_2UE=sum(C_2UE_real)/R
CSI_2UE=sum(CSI_2UE_real)/R