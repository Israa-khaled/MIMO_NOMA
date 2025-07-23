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
MU=20;% the users belongs to the same side of BS, i.e., 0<theta<180.
r=1;R=1;beta_seuil=0.5;TXPower_mw=Pe*1e3;TXPower=10*log10(TXPower_mw);NP=noise_power;
gamma_FPA1=0.4;
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
%2-UE clustering:
[Clusters_array_2UE,nb_2UEinCluster_2UE,nb_clusters_2UE]=get_UserClustering_2UE(beta_n_seuil,MU);
%% The beam associated to each cluster
% Associate a beam for each cluster by using DBS
% 2-UE NOMA-DBS
[WP_DBSn_noma_2UE,W_DBSn_noma_2UE]= getWdbs_noma_multiUE1(Clusters_array_2UE,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_2UE,nb_clusters_2UE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
aux_mat_allUEs_noma_2UE= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma_2UE,nb_clusters_2UE);
%Partial-CSI
aux_mat_allUEs_noma_2UE_partialCSI= getAuxMatrix_noma_partialCSI(MU,a_BS,WP_DBSn_noma_2UE,nb_clusters_2UE);
%% Zeta-based User Ordering
% FCSI-FPA with gamma1
[P_noma_2UE_FPA1, rate_user_NOMA_FCT_2UE_FPA1]=  getGamma_FPA_2UE(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,gamma_FPA1,NP);
aux_mat_allUEs_2UE_FPA1= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA1,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA1= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA1,aux_mat_allUEs_2UE_FPA1,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA1(r)=sum(rate_user_noma_2UE_FPA1);
% ADI-FPA with gamma1
[P_noma_2UE_FPA1_ADI, rate_user_NOMA_FCT_2UE_FPA1_ADI]=  getGamma_FPA_2UE_ADI(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,aux_mat_allUEs_noma_2UE_partialCSI,gamma_FPA1,NP);
aux_mat_allUEs_2UE_FPA1_ADI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA1_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA1_ADI= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA1_ADI,aux_mat_allUEs_2UE_FPA1_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA1_ADI(r)=sum(rate_user_noma_2UE_FPA1_ADI);
%% Distance-based User Ordering 
% FPA with gamma1
[P_noma_2UE_FPA1_dist, rate_user_NOMA_FCT_2UE_FPA1_dist]=  getGamma_FPA_2UE_Distance(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,gamma_FPA1,NP,TR_dist);
aux_mat_allUEs_2UE_FPA1= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA1_dist,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA1_dist= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA1_dist,aux_mat_allUEs_2UE_FPA1,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA1_dist(r)=sum(rate_user_noma_2UE_FPA1_dist);
%% HW-based User Ordering
% FPA with gamma1
[P_noma_2UE_FPA1_ADI_HW, rate_user_NOMA_FCT_2UE_FPA1_ADI_HW]=  getGamma_FPA_2UE_HW(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,gamma_FPA1,NP,HW_DBS);
aux_mat_allUEs_2UE_FPA1_ADI_HW= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA1_ADI_HW,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA1_ADI_HW= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA1_ADI_HW,aux_mat_allUEs_2UE_FPA1_ADI_HW,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA1_ADI_HW(r)=sum(rate_user_noma_2UE_FPA1_ADI_HW);
%%
r=r+1;
end
%% average SUMRATE
sumrate_DBS=sum(sumrate_DBS1)/R;
sumrate_CB=sum(sumrate_CB1)/R;
sumrate_ZF=sum(sumrate_ZF1)/R;

% full csi
sumrate_dbsnoma_2UE_FPA1=sum(sumrate_dbsnoma_2UE_real_FPA1)/R;
% partial csi 
sumrate_dbsnoma_2UE_FPA1_ADI=sum(sumrate_dbsnoma_2UE_real_FPA1_ADI)/R;

%
sumrate_dbsnoma_2UE_FPA1_ADI_HW=sum(sumrate_dbsnoma_2UE_real_FPA1_ADI_HW)/R;
sumrate_dbsnoma_2UE_FPA1_dist=sum(sumrate_dbsnoma_2UE_real_FPA1_dist)/R;

