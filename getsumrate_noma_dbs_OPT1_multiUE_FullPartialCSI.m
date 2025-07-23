function [sumrate_DBS,sumrate_CB,sumrate_ZF,sumrate_dbsnoma_multiUE,sumrate_dbsnoma_2UE]= getsumrate_noma_dbs_OPT1_multiUE_FullPartialCSI(MU,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe)
addpath TempFct2UE mUE
rmpath nyusimv1.6.1
if scn == 'mono'
  addpath matlab_code_NYUSIM_monotraject
  rmpath matlab_code_NYUSIM_multitraject  
elseif scn == 'mult'
  rmpath matlab_code_NYUSIM_monotraject
  addpath matlab_code_NYUSIM_multitraject
end
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
load inputParameters 
r=1;TXPower_mw=Pe*1e3;TXPower=10*log10(TXPower_mw);NP=noise_power;
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
% multi-UE NOMA-DBS
[WP_DBSn_noma_multiUE,W_DBSn_noma_multiUE]= getWdbs_noma_multiUE1(Clusters_array_multiUE,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_multiUE,nb_clusters_multiUE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
aux_mat_allUEs_noma_multiUE= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma_multiUE,nb_clusters_multiUE);
% 2-UE NOMA-DBS
[WP_DBSn_noma_2UE,W_DBSn_noma_2UE]= getWdbs_noma_multiUE1(Clusters_array_2UE,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_2UE,nb_clusters_2UE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
aux_mat_allUEs_noma_2UE= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma_2UE,nb_clusters_2UE);
%%  Power allocation method
% multi-UE NOMA-DBS
[P_noma_multiUE, rate_user_NOMA_FCT_multiUE]= getGamma_OptPA_multiUE1(Clusters_array_multiUE,MU,nb_2UEinCluster_multiUE,nb_clusters_multiUE,HW_DBS,aux_mat_allUEs_noma_multiUE,Pmin,Pe,NP);
% 2-UE NOMA-DBS
[P_noma_2UE, rate_user_NOMA_FCT_2UE]=  getGamma_OptPA_multiUE1(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,HW_DBS,aux_mat_allUEs_noma_2UE,Pmin,Pe,NP);
%% Sum-rate evaluation
%AUXILIAR MATRIX to calculate the SINR of each user 
% multi-UE NOMA-DBS
aux_mat_allUEs_multiUE= getAuxMatrix_multiUE1(Clusters_array_multiUE,MU,H_total1,WP_DBSn_noma_multiUE,P_noma_multiUE,nb_clusters_multiUE,nb_2UEinCluster_multiUE);
rate_user_noma_multiUE= getRate_DBS_NOMA_multiUE(Clusters_array_multiUE,rate_user_NOMA_FCT_multiUE,aux_mat_allUEs_multiUE,nb_clusters_multiUE,nb_2UEinCluster_multiUE,HW_DBS,NP);
sumrate_dbsnoma_multiUE_real(r)=sum(rate_user_noma_multiUE);
% 2-UE NOMA-DBS
aux_mat_allUEs_2UE= getAuxMatrix_multiUE1(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2UE,P_noma_2UE,nb_clusters_2UE,nb_2UEinCluster_2UE);
rate_user_noma_2UE= getRate_DBS_NOMA_multiUE(Clusters_array_2UE,rate_user_NOMA_FCT_2UE,aux_mat_allUEs_2UE,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2UE_real(r)=sum(rate_user_noma_2UE);

%%
r=r+1;
end
%% average SUMRATE
sumrate_DBS=sum(sumrate_DBS1)/R;
sumrate_CB=sum(sumrate_CB1)/R;
sumrate_ZF=sum(sumrate_ZF1)/R;
sumrate_dbsnoma_multiUE=sum(sumrate_dbsnoma_multiUE_real)/R;
sumrate_dbsnoma_2UE=sum(sumrate_dbsnoma_2UE_real)/R;