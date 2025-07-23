clear all;close all;clc;
%%israa
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180. 
f=28;sceType = 'UMi';Pe=1;
Nt=64;TxArrayType='URA';Wt=32;scn='mono';
MU=7;% the users belongs to the same side of BS, i.e., 0<theta<180.
r=1;R=1;beta_seuil=0.5;
rmpath nyusimv1.6.1
if scn == 'mono'
  addpath matlab_code_NYUSIM_monotraject
  rmpath matlabcode_article1  
elseif scn == 'mult'
  rmpath matlab_code_NYUSIM_monotraject
  addpath matlabcode_article1 
end
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
load inputParameters 
RFBW=20;Tech=1000/(RFBW);%en ns 
wl=3e8/(f*1e9);
r=1;
addpath TempFct
if scn == 'mono'
  addpath matlab_code_NYUSIM_monotraject
  rmpath matlabcode_article1  
elseif scn == 'mult'
  rmpath matlab_code_NYUSIM_monotraject
  addpath matlabcode_article1 
end
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
load inputParameters 
r=1;
while (r<=R)
[CIR_MIMO_Struct,powerSpectrum,H_total1,H,AOD_LobePowerSpectrum,powerspectrum_struct,Pr_dBm,TR_dist] = getH_MIMO_n_1(f,RFBW,sceType,envType,TXPower,Nt,Nr,MU,dmin,dmax,h_BS,TxArrayType,Wt,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX);  
n_BS = [1:1:Nt]';
j_complex = sqrt(-1);
% DBS preocding matrix
[azimuthAOD_USER,eleAOD_USER,a_BS,W_DBSn]= getWdbs(powerspectrum_struct,Nt,Wt,TxArrayType,MU,dTxAnt);  
HW_DBS=H_total1*W_DBSn;
[SINR_user_DBS,SIR_user_DBS,rate_user_DBS] =getsumrateULAURA_Channelisgiven(a_BS,H_total1,MU,noise_power);
sumrate_DBS1(r)=sum(rate_user_DBS);

%% Calculate beta matrix (MU-1,MU)
% group all the user that has beta >seuil 
[beta_n_temp]= betaMatrix_fctbeta0(a_BS,beta_seuil,MU,Nt);

%% Group users two-by-two 
[nb_2UEinCluster,nb_clusters,Clusters_array]= getClustersArray(beta_n_temp,MU);
%% The beam associated to each cluster
% Associate a beam beam for each cluster by using DBS
% DBS with NOMA
[WP_DBSn_noma,W_DBSn_noma]= getWdbs_noma(Clusters_array,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster,nb_clusters,TxArrayType,Nt,Wt,dTxAnt,Pe,MU);

%power allocated to each user 
[P_noma]= getGamma_classicalPA(Clusters_array,MU,nb_2UEinCluster,HW_DBS);
%AUXILIAR MATRIX to calculate the SINR of each user 
aux_mat_allUEs_UnifPA= getAuxMatrix(Clusters_array,MU,H_total1,W_DBSn_noma,P_noma,nb_clusters,nb_2UEinCluster);% For uniform PA
aux_mat_allUEs_NUnifPA= getAuxMatrix(Clusters_array,MU,H_total1,WP_DBSn_noma,P_noma,nb_clusters,nb_2UEinCluster);% for no uniform PA
% rate of each user
NP=noise_power;
rate_user_noma_uPA= getRate_DBS_NOMA(Clusters_array,aux_mat_allUEs_UnifPA,nb_clusters,nb_2UEinCluster,HW_DBS,NP);
sumrate_dbsnoma1_uPA(r)=sum(rate_user_noma_uPA);
rate_user_noma_nuPA= getRate_DBS_NOMA(Clusters_array,aux_mat_allUEs_NUnifPA,nb_clusters,nb_2UEinCluster,HW_DBS,NP);
sumrate_dbsnoma1_nuPA(r)=sum(rate_user_noma_nuPA);

r=r+1;
end
%% SUMRATE for all realizations
sumrate_DBS=sum(sumrate_DBS1)/R;
sumrate_dbsnoma_uPA=sum(sumrate_dbsnoma1_uPA)/R;
sumrate_dbsnoma_nuPA=sum(sumrate_dbsnoma1_nuPA)/R;

