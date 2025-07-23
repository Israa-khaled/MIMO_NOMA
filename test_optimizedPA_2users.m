clear all; clc;
rmpath matlab_code_NYUSIM_monotraject
rmpath matlabcode_article1 
scn='mono'; MU=5;Nt=64;Wt=1;sceType='UMi';TxArrayType='ULA';R=1;f=28;beta_seuil=0.4;TXPower=1;
addpath nyusimv1.6.1 
% N.B. the users belongs to the same side of BS, i.e., 0<theta<180.
load inputParameters 
RFBW=20;Tech=1000/(RFBW);%en ns 
wl=3e8/(f*1e9);
r=1;
% 2 USERS 
theta_user1=30;phi_user1=0;TRDistance1=10;
theta_user2=32;phi_user2=0;TRDistance2=100;
[H_ensemble1,powerspectrum_struct1,TRdist1] = w_getH_MIMO_userdirection_1UE(f,RFBW,sceType,envType,dmin,dmax,TXPower,h_BS,Nt,Nr,1,TxArrayType,RxArrayType,Wt,Wr,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_user1,phi_user1,TRDistance1);
[H_ensemble2,powerspectrum_struct2,TRdist2] = w_getH_MIMO_userdirection_1UE(f,RFBW,sceType,envType,dmin,dmax,TXPower,h_BS,Nt,Nr,1,TxArrayType,RxArrayType,Wt,Wr,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_user2,phi_user2,TRDistance2);
%% DBS 
%User 1
powerSpectrum=powerspectrum_struct1.(['powerspectrum_',num2str(1)]);
aziAOD1 = powerSpectrum(1,4) 
eleAOD = powerSpectrum(1,5);
n_BS = [1:1:Nt]';
j_complex = sqrt(-1);
for n_bs = 1:length(n_BS)
  a_BS(n_bs,1) = exp(j_complex.*(n_BS(n_bs)-1).*2.*pi.*dTxAnt.*cos(aziAOD1*pi/180));    
end
%User 2
powerSpectrum=powerspectrum_struct2.(['powerspectrum_',num2str(1)]);
aziAOD2 = powerSpectrum(1,4) 
eleAOD = powerSpectrum(1,5);
n_BS = [1:1:Nt]';
j_complex = sqrt(-1);
for n_bs = 1:length(n_BS)
  a_BS(n_bs,2) = exp(j_complex.*(n_BS(n_bs)-1).*2.*pi.*dTxAnt.*cos(aziAOD2*pi/180));    
end
H_TOT(1,:)=H_ensemble1;H_TOT(2,:)=H_ensemble2;
[SINR_user_DBS,SIR_user_DBS,rate_user_DBS] = getsumrateULAURA_Channelisgiven(a_BS,H_TOT,2,noise_power);

%% NOMA
beta=abs(a_BS(:,1)'*a_BS(:,2))/Nt;
% beam cluster for NOMA-DBS 
aziAOD=(aziAOD1+aziAOD2)/2
for n_bs = 1:length(n_BS)
    a_BS_ULA_Noma(n_bs,1) = exp(j_complex.*(n_BS(n_bs)-1).*2.*pi.*dTxAnt.*cos(aziAOD*pi/180));    
end
W_DBS=a_BS;
W_DBSn=W_DBS/sqrt(trace(W_DBS*W_DBS'));
% power allocation for NOMA-DBS
% CASE 1: typical PA
HW_DBS=H_TOT*W_DBSn;
h1w1=abs(HW_DBS(1,1))^2;
h2w2=abs(HW_DBS(2,2))^2;
P_noma(1)= h2w2/(h1w1+h2w2);
P_noma(2)= h1w1/(h1w1+h2w2);
% CASE 2: Optimized PA
h1w=H_TOT(1,:)*W_DBSn;
h2w=H_TOT(2,:)*W_DBSn;
w1w=abs(a_BS(:,1)'*a_BS_ULA_Noma)/64;
w2w=abs(a_BS(:,2)'*a_BS_ULA_Noma)/64;
