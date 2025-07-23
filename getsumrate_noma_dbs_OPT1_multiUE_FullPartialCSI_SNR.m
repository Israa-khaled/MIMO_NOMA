function [sumrate_DBS,sumrate_CB,sumrate_ZF,sumrate_dbsnoma_multiUE_fullCSI,sumrate_dbsnoma_2UE_fullCSI,sumrate_dbsnoma_multiUE_partialCSI,sumrate_dbsnoma_2UE_partialCSI]= getsumrate_noma_dbs_OPT1_multiUE_FullPartialCSI_SNR(MU,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe,NP)
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
r=1;TXPower_mw=Pe*1e3;TXPower=10*log10(TXPower_mw);
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
NbUEsinCluster_array=zeros(MU,R);
while (r<=R)
[CIR_MIMO_Struct,powerSpectrum,H_total1,H,AOD_LobePowerSpectrum,powerspectrum_struct,Pr_dBm,TR_dist] = getH_MIMO_n_1(f,RFBW,sceType,envType,TXPower,Nt,Nr,MU,dmin,dmax,h_BS,TxArrayType,Wt,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX);  
% DBS preocding matrix
[azimuthAOD_USER,eleAOD_USER,a_BS,W_DBSn]= getWdbs(powerspectrum_struct,Nt,Wt,TxArrayType,MU,dTxAnt);  
HW_DBS=H_total1*W_DBSn;
[SINR_user_DBS,SINR_user_CB,SINR_user_ZF,rate_user_DBS,rate_user_CB,rate_user_ZF] =getsumrateULAURA_Channelisgiven(a_BS,H_total1,MU,NP);
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
sum_UEsinmultiUEcluster=0;
for index_c=1:nb_2UEinCluster_multiUE
    sum_UEsinmultiUEcluster(index_c)=length(Clusters_array_multiUE{index_c});
end
Average_sum_UEsinmultiUEcluster(r)=prctile(sum_UEsinmultiUEcluster,50);
r=r+1;

end
%% average SUMRATE
sumrate_DBS=sum(sumrate_DBS1)/R;
sumrate_CB=sum(sumrate_CB1)/R;
sumrate_ZF=sum(sumrate_ZF1)/R;
% full csi
sumrate_dbsnoma_multiUE_fullCSI=sum(sumrate_dbsnoma_multiUE_real_fullCSI)/R;
sumrate_dbsnoma_2UE_fullCSI=sum(sumrate_dbsnoma_2UE_real_fullCSI)/R;
% partial csi 
sumrate_dbsnoma_multiUE_partialCSI=sum(sumrate_dbsnoma_multiUE_real_partialCSI)/R;
sumrate_dbsnoma_2UE_partialCSI=sum(sumrate_dbsnoma_2UE_real_partialCSI)/R;
Average_sum_UEsinmultiUEcluster_allReal=sum(Average_sum_UEsinmultiUEcluster)/R;