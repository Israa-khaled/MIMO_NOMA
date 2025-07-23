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
gamma_FPA1=0.1;gamma_FPA2=0.2;gamma_FPA3=0.5;
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
[Clusters_array_multiUE,nb_2UEinCluster_multiUE,nb_clusters_multiUE]=get_UserClustering_multiUE(MU,Grp2b2Interf,beta_Grp2b2Interf);
%2-UE clustering:
[Clusters_array_2UE,nb_2UEinCluster_2UE,nb_clusters_2UE]=get_UserClustering_2UE(beta_n_seuil,MU);
%% The beam associated to each cluster
% Associate a beam for each cluster by using DBS
%% multi-UE NOMA-DBS
[WP_DBSn_noma_multiUE,W_DBSn_noma_multiUE]= getWdbs_noma_multiUE1(Clusters_array_multiUE,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_multiUE,nb_clusters_multiUE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
%Full-CSI
aux_mat_allUEs_noma_multiUE_fullCSI= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma_multiUE,nb_clusters_multiUE);
%Partial-CSI
aux_mat_allUEs_noma_multiUE_partialCSI= getAuxMatrix_noma_partialCSI(MU,a_BS,WP_DBSn_noma_multiUE,nb_clusters_multiUE);
%% 2-UE NOMA-DBS
[WP_DBSn_noma_2UE,W_DBSn_noma_2UE]= getWdbs_noma_multiUE1(Clusters_array_2UE,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_2UE,nb_clusters_2UE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
% FULL-CSI 
aux_mat_allUEs_noma_2UE_fullCSI= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma_2UE,nb_clusters_2UE);
%Partial-CSI
aux_mat_allUEs_noma_2UE_partialCSI= getAuxMatrix_noma_partialCSI(MU,a_BS,WP_DBSn_noma_2UE,nb_clusters_2UE);
%%  Power allocation method
% multi-UE NOMA-DBS||Full-CSI
[P_noma_multiUE_fullCSI, rate_user_NOMA_FCT_multiUE_fullCSI]=getGamma_OptPA_multiUE1(Clusters_array_multiUE,MU,nb_2UEinCluster_multiUE,nb_clusters_multiUE,HW_DBS,aux_mat_allUEs_noma_multiUE_fullCSI,Pmin,Pe,NP);
% multi-UE NOMA-DBS||Partial-CSI
[P_noma_multiUE_partialCSI, rate_user_NOMA_FCT_multiUE_partialCSI]=getGamma_OptPA_multiUE1_partialCSI(Clusters_array_multiUE,MU,nb_2UEinCluster_multiUE,nb_clusters_multiUE,HW_DBS,aux_mat_allUEs_noma_multiUE_fullCSI,aux_mat_allUEs_noma_multiUE_partialCSI,Pmin,Pe,NP);
% 2-UE NOMA-DBS||Full-CSI
[P_noma_2UE_fullCSI, rate_user_NOMA_FCT_2UE_fullCSI]=getGamma_OptPA_multiUE1(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,HW_DBS,aux_mat_allUEs_noma_2UE_fullCSI,Pmin,Pe,NP);
% 2-UE NOMA-DBS||Partial-CSI
[P_noma_2UE_partialCSI, rate_user_NOMA_FCT_2UE_partialCSI]=getGamma_OptPA_multiUE1_partialCSI(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,HW_DBS,aux_mat_allUEs_noma_2UE_fullCSI,aux_mat_allUEs_noma_2UE_partialCSI,Pmin,Pe,NP);
%% Sum-rate evaluation
%AUXILIAR MATRIX to calculate the SINR of each user 
% multi-UE NOMA-DBS||Full-CSI
aux_mat_allUEs_multiUE_fullCSI= getAuxMatrix_multiUE1(Clusters_array_multiUE,MU,H_total1,WP_DBSn_noma_multiUE,P_noma_multiUE_fullCSI,nb_clusters_multiUE,nb_2UEinCluster_multiUE);
rate_user_noma_multiUE_fullCSI= getRate_DBS_NOMA_multiUE(Clusters_array_multiUE,rate_user_NOMA_FCT_multiUE_fullCSI,aux_mat_allUEs_multiUE_fullCSI,nb_clusters_multiUE,nb_2UEinCluster_multiUE,HW_DBS,NP);
sumrate_dbsnoma_multiUE_real_fullCSI(r)=sum(rate_user_noma_multiUE_fullCSI);
% multi-UE NOMA-DBS||Partial-CSI
aux_mat_allUEs_multiUE_partialCSI= getAuxMatrix_multiUE1(Clusters_array_multiUE,MU,H_total1,WP_DBSn_noma_multiUE,P_noma_multiUE_partialCSI,nb_clusters_multiUE,nb_2UEinCluster_multiUE);
rate_user_noma_multiUE_partialCSI= getRate_DBS_NOMA_multiUE(Clusters_array_multiUE,rate_user_NOMA_FCT_multiUE_partialCSI,aux_mat_allUEs_multiUE_partialCSI,nb_clusters_multiUE,nb_2UEinCluster_multiUE,HW_DBS,NP);
sumrate_dbsnoma_multiUE_real_partialCSI(r)=sum(rate_user_noma_multiUE_partialCSI);
% 2-UE NOMA-DBS||Full-CSI
aux_mat_allUEs_2UE_fullCSI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_fullCSI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_fullCSI= getRate_DBS_NOMA_multiUE(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_fullCSI,aux_mat_allUEs_2UE_fullCSI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_fullCSI(r)=sum(rate_user_noma_2UE_fullCSI);
% 2-UE NOMA-DBS||Partial-CSI
aux_mat_allUEs_2UE_partialCSI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_partialCSI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_partialCSI= getRate_DBS_NOMA_multiUE(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_partialCSI,aux_mat_allUEs_2UE_partialCSI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_partialCSI(r)=sum(rate_user_noma_2UE_partialCSI);
%%
%% The beam associated to each cluster
% Associate a beam for each cluster by using DBS
% 2-UE NOMA-DBS
[WP_DBSn_noma_2UE,W_DBSn_noma_2UE]= getWdbs_noma_multiUE1(Clusters_array_2UE,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_2UE,nb_clusters_2UE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
aux_mat_allUEs_noma_2UE= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma_2UE,nb_clusters_2UE);
%Partial-CSI
aux_mat_allUEs_noma_2UE_partialCSI= getAuxMatrix_noma_partialCSI(MU,a_BS,WP_DBSn_noma_2UE,nb_clusters_2UE);
% Power allocation method: FPA 
% FPA with gamma1=0.1
[P_noma_2UE_FPA1, rate_user_NOMA_FCT_2UE_FPA1]=  getGamma_FPA_2UE(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,gamma_FPA1,NP);
aux_mat_allUEs_2UE_FPA1= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA1,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA1= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA1,aux_mat_allUEs_2UE_FPA1,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA1(r)=sum(rate_user_noma_2UE_FPA1);

% FPA with gamma1=0.2
[P_noma_2UE_FPA2, rate_user_NOMA_FCT_2UE_FPA2]=  getGamma_FPA_2UE(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,gamma_FPA2,NP);
aux_mat_allUEs_2UE_FPA2= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA2,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA2= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA2,aux_mat_allUEs_2UE_FPA2,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA2(r)=sum(rate_user_noma_2UE_FPA2);

% FPA with gamma1=0.3
[P_noma_2UE_FPA3, rate_user_NOMA_FCT_2UE_FPA3]=  getGamma_FPA_2UE(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,gamma_FPA3,NP);
aux_mat_allUEs_2UE_FPA3= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA3,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA3= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA3,aux_mat_allUEs_2UE_FPA3,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA3(r)=sum(rate_user_noma_2UE_FPA3);

% FTPA with tau1=0.2
tau=0.2;
[P_noma_2UE_FTPA, rate_user_NOMA_FCT_2UE_FTPA]=  getGamma_FTPA_2UE(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,HW_DBS,NP,tau);
aux_mat_allUEs_2UE_FTPA= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FTPA,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FTPA= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FTPA,aux_mat_allUEs_2UE_FTPA,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FTPA(r)=sum(rate_user_noma_2UE_FTPA);

% FTPA with tau1=0.5
tau=0.5;
[P_noma_2UE_FTPA, rate_user_NOMA_FCT_2UE_FTPA]=  getGamma_FTPA_2UE(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,HW_DBS,NP,tau);
aux_mat_allUEs_2UE_FTPA= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FTPA,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FTPA= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FTPA,aux_mat_allUEs_2UE_FTPA,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FTPA1(r)=sum(rate_user_noma_2UE_FTPA);

% Power allocation method: FPA 
% FPA with gamma1=0.1
[P_noma_2UE_FPA1_ADI, rate_user_NOMA_FCT_2UE_FPA1_ADI]=  getGamma_FPA_2UE_ADI(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,aux_mat_allUEs_noma_2UE_partialCSI,gamma_FPA1,NP);
aux_mat_allUEs_2UE_FPA1_ADI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA1_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA1_ADI= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA1_ADI,aux_mat_allUEs_2UE_FPA1_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA1_ADI(r)=sum(rate_user_noma_2UE_FPA1_ADI);

% FPA with gamma1=0.2
[P_noma_2UE_FPA2_ADI, rate_user_NOMA_FCT_2UE_FPA2_ADI]=  getGamma_FPA_2UE_ADI(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,aux_mat_allUEs_noma_2UE_partialCSI,gamma_FPA2,NP);
aux_mat_allUEs_2UE_FPA2_ADI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA2_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA2_ADI= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA2_ADI,aux_mat_allUEs_2UE_FPA2_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA2_ADI(r)=sum(rate_user_noma_2UE_FPA2_ADI);

% FPA with gamma1=0.3
[P_noma_2UE_FPA3_ADI, rate_user_NOMA_FCT_2UE_FPA3_ADI]=  getGamma_FPA_2UE_ADI(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,aux_mat_allUEs_noma_2UE_partialCSI,gamma_FPA3,NP);
aux_mat_allUEs_2UE_FPA3_ADI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FPA3_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FPA3_ADI= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FPA3_ADI,aux_mat_allUEs_2UE_FPA3_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FPA3_ADI(r)=sum(rate_user_noma_2UE_FPA3_ADI);

% FTPA 
tau=0.000002;
[P_noma_2UE_FTPA_ADI, rate_user_NOMA_FCT_2UE_FTPA_ADI]=  getGamma_FTPA_2UE_ADI(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,aux_mat_allUEs_noma_2UE_partialCSI,HW_DBS,NP,tau);
aux_mat_allUEs_2UE_FTPA_ADI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FTPA_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FTPA_ADI= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FTPA_ADI,aux_mat_allUEs_2UE_FTPA_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FTPA_ADI(r)=sum(rate_user_noma_2UE_FTPA_ADI);

% FTPA 
tau=0.5;
[P_noma_2UE_FTPA_ADI, rate_user_NOMA_FCT_2UE_FTPA_ADI]=  getGamma_FTPA_2UE_ADI(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,aux_mat_allUEs_noma_2UE,aux_mat_allUEs_noma_2UE_partialCSI,HW_DBS,NP,tau);
aux_mat_allUEs_2UE_FTPA_ADI= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE_FTPA_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE_FTPA_ADI= getRate_DBS_NOMA_multiUE_v1(Clusters_array_2UE,rate_user_NOMA_FCT_2UE_FTPA_ADI,aux_mat_allUEs_2UE_FTPA_ADI,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real_FTPA_ADI1(r)=sum(rate_user_noma_2UE_FTPA_ADI);
%%
r=r+1;
end
%% average SUMRATE
sumrate_DBS=sum(sumrate_DBS1)/R;
sumrate_CB=sum(sumrate_CB1)/R;
sumrate_ZF=sum(sumrate_ZF1)/R;

% full csi
sumrate_dbsnoma_multiUE_fullCSI=sum(sumrate_dbsnoma_multiUE_real_fullCSI)/R;
sumrate_dbsnoma_2UE_fullCSI=sum(sumrate_dbsnoma_2UE_real_fullCSI)/R;

sumrate_dbsnoma_2UE_FPA1=sum(sumrate_dbsnoma_2UE_real_FPA1)/R;
sumrate_dbsnoma_2UE_FPA2=sum(sumrate_dbsnoma_2UE_real_FPA2)/R;
sumrate_dbsnoma_2UE_FPA3=sum(sumrate_dbsnoma_2UE_real_FPA3)/R;
sumrate_dbsnoma_2UE_FTPA=sum(sumrate_dbsnoma_2UE_real_FTPA)/R;
sumrate_dbsnoma_2UE_FTPA1=sum(sumrate_dbsnoma_2UE_real_FTPA1)/R;

% partial csi 
sumrate_dbsnoma_multiUE_partialCSI=sum(sumrate_dbsnoma_multiUE_real_partialCSI)/R;
sumrate_dbsnoma_2UE_partialCSI=sum(sumrate_dbsnoma_2UE_real_partialCSI)/R;
sumrate_dbsnoma_2UE_FPA1_ADI=sum(sumrate_dbsnoma_2UE_real_FPA1_ADI)/R;
sumrate_dbsnoma_2UE_FPA2_ADI=sum(sumrate_dbsnoma_2UE_real_FPA2_ADI)/R;
sumrate_dbsnoma_2UE_FPA3_ADI=sum(sumrate_dbsnoma_2UE_real_FPA3_ADI)/R;
sumrate_dbsnoma_2UE_FTPA_ADI=sum(sumrate_dbsnoma_2UE_real_FTPA_ADI)/R;
sumrate_dbsnoma_2UE_FTPA_ADI1=sum(sumrate_dbsnoma_2UE_real_FTPA_ADI1)/R;

