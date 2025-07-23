clear all; clc;
rmpath matlab_code_NYUSIM_monotraject
rmpath matlabcode_article1 
f=28;sceType = 'UMi';Pe=1;Pmin=1e-3;
Nt=32;TxArrayType='ULA';Wt=1;scn='mono';
MU=4;% the users belongs to the same side of BS, i.e., 0<theta<180.
r=1;R=1;beta_seuil=0.5;TXPower=Pe;
addpath nyusimv1.6.1
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
load inputParameters 
RFBW=20;Tech=1000/(RFBW);%en ns 
wl=3e8/(f*1e9);
r=1;
% 2 USERS 
theta_user1=30.5;phi_user1=0;TRDistance1=10;
theta_user2=30;phi_user2=0;TRDistance2=100;
theta_user3=31;phi_user3=0;TRDistance3=30;
theta_user4=32;phi_user4=0;TRDistance4=60;
[H_ensemble1,powerspectrum_struct1,TRdist1] = w_getH_MIMO_userdirection_1UE(f,RFBW,sceType,envType,dmin,dmax,TXPower,h_BS,Nt,Nr,1,TxArrayType,RxArrayType,Wt,Wr,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_user1,phi_user1,TRDistance1);
[H_ensemble2,powerspectrum_struct2,TRdist2] = w_getH_MIMO_userdirection_1UE(f,RFBW,sceType,envType,dmin,dmax,TXPower,h_BS,Nt,Nr,1,TxArrayType,RxArrayType,Wt,Wr,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_user2,phi_user2,TRDistance2);
[H_ensemble3,powerspectrum_struct3,TRdist3] = w_getH_MIMO_userdirection_1UE(f,RFBW,sceType,envType,dmin,dmax,TXPower,h_BS,Nt,Nr,1,TxArrayType,RxArrayType,Wt,Wr,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_user3,phi_user1,TRDistance3);
[H_ensemble4,powerspectrum_struct4,TRdist4] = w_getH_MIMO_userdirection_1UE(f,RFBW,sceType,envType,dmin,dmax,TXPower,h_BS,Nt,Nr,1,TxArrayType,RxArrayType,Wt,Wr,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_user4,phi_user2,TRDistance4);
n_BS = [1:1:Nt]'; H_TOT(1,:)=H_ensemble1;H_TOT(2,:)=H_ensemble2;H_TOT(3,:)=H_ensemble3;H_TOT(4,:)=H_ensemble4;
j_complex = sqrt(-1);H_total1=H_TOT;
powerSpectrum=powerspectrum_struct1.(['powerspectrum_',num2str(1)]);
powerspectrum_struct.(['powerspectrum_',num2str(1)])=powerSpectrum;
powerSpectrum=powerspectrum_struct2.(['powerspectrum_',num2str(1)]);
powerspectrum_struct.(['powerspectrum_',num2str(2)])=powerSpectrum;
powerSpectrum=powerspectrum_struct3.(['powerspectrum_',num2str(1)]);
powerspectrum_struct.(['powerspectrum_',num2str(3)])=powerSpectrum;
powerSpectrum=powerspectrum_struct4.(['powerspectrum_',num2str(1)]);
powerspectrum_struct.(['powerspectrum_',num2str(4)])=powerSpectrum;
% powerspectrum_struct.(['powerspectrum_',num2str(2)])=powerSpectrum_struct2.(['powerspectrum_',num2str(1)]);
% powerspectrum_struct.(['powerspectrum_',num2str(3)])=powerSpectrum_struct3.(['powerspectrum_',num2str(1)]);
% powerspectrum_struct.(['powerspectrum_',num2str(4)])=powerSpectrum_struct4.(['powerspectrum_',num2str(1)]);
% DBS preocding matrix
[azimuthAOD_USER,eleAOD_USER,a_BS,W_DBSn]= getWdbs(powerspectrum_struct,Nt,Wt,TxArrayType,MU,dTxAnt);  
HW_DBS=H_total1*W_DBSn;
[SINR_user_DBS,SINR_user_CB,SINR_user_ZF,rate_user_DBS,rate_user_CB,rate_user_ZF] =getsumrateULAURA_Channelisgiven(a_BS,H_total1,MU,noise_power);
sumrate_DBS1=sum(rate_user_DBS)
%% User Clustering Algorithm
% Select the set of groups of 2-UEs who have a spatial interference greater than  \beta_0 
[Grp2b2Interf,beta_n_allUser,beta_n_seuil,beta_Grp2b2Interf]= UserbetaMatrix_beta0(a_BS,beta_seuil,MU,Nt);
%Multi-UE clustering:
[Clusters_array,nb_2UEinCluster,nb_clusters]=get_UserClustering_multiUE(MU,Grp2b2Interf,beta_Grp2b2Interf);
%2-UE clustering:
%% Group users two-by-two 
NP=noise_power;
[nb_2UEinCluster_2UE,nb_clusters_2UE,Clusters_array_2UE]= getClustersArray(beta_n_seuil,MU);
% [sumrate_dbsnoma1_uPA,sumrate_dbsnoma1_nuPA] =getAll(H_total1,Clusters_array_2UE,HW_DBS,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_2UE,nb_clusters_2UE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU,noise_power)
%% NOMA-DBS 2UE 
% The beam associated to each cluster
% Associate a beam beam for each cluster by using DBS
% DBS with NOMA
[WP_DBSn_noma_2ue,W_DBSn_noma_2ue]= getWdbs_noma(Clusters_array_2UE,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster_2UE,nb_clusters_2UE,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
aux_mat_allUEs_noma_2ue= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma_2ue,nb_clusters_2UE);
%power allocated to each user 
NP=noise_power;
[P_noma_SINR, rate_user_NOMA_FCT_SINR,nbuserActive_NOMA_SINR]= getGamma_OptPA_withoutConstraints(Clusters_array_2UE,MU,nb_2UEinCluster_2UE,nb_clusters_2UE,HW_DBS,aux_mat_allUEs_noma_2ue,Pmin,Pe,NP);
%AUXILIAR MATRIX to calculate the SINR of each user 
aux_mat_allUEs_2UE= getAuxMatrix(Clusters_array_2UE,MU,H_total1,WP_DBSn_noma_2ue,P_noma_SINR,nb_clusters_2UE,nb_2UEinCluster_2UE);
% rate of each user
rate_user_noma_2UE =getRate_DBS_NOMA(Clusters_array_2UE,aux_mat_allUEs_2UE,nb_clusters_2UE,nb_2UEinCluster_2UE,HW_DBS,NP);
sumrate_dbsnoma_2ue=sum(rate_user_noma_2UE)
%%NOMA-DBS multi-UE
[WP_DBSn_noma,W_DBSn_noma]= getWdbs_noma_multiUE1(Clusters_array,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster,nb_clusters,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);
aux_mat_allUEs_noma= getAuxMatrix_noma(MU,H_total1,WP_DBSn_noma,nb_clusters);
NP=noise_power;
[P_noma, rate_user_NOMA_FCT]= getGamma_OptPA_multiUE(Clusters_array,MU,nb_2UEinCluster,nb_clusters,HW_DBS,aux_mat_allUEs_noma,Pmin,Pe,NP);
aux_mat_allUEs= getAuxMatrix_multiUE(Clusters_array,MU,H_total1,W_DBSn_noma,P_noma,nb_clusters,nb_2UEinCluster);
rate_user_noma_multiUE= getRate_DBS_NOMA_multiUE(Clusters_array,rate_user_NOMA_FCT,aux_mat_allUEs,nb_clusters,nb_2UEinCluster,HW_DBS,NP);
sumrate_dbsnoma_multiUE=sum(rate_user_noma_multiUE)
clear rate_user_NOMA_FCT_SINR  rate_user_NOMA_FCT