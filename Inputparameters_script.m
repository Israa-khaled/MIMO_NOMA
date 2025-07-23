clear all,clc,
%  -donnees ---------------------------------------------------------------
%f = 28;% en GHz
RFBW = 20;% en MHz
% Tech=1000/(RFBW);%en ns 
%TXPower = 30;% en dBm 
%Pe=1;
% Transmit antenna spacing in wavelengths (0.1-100)
dTxAnt = 0.5; 
% Receive antenna spacing in wavelengths (0.1-100)
dRxAnt = 0.5;
% Transmit antenna azimuth half-power beamwidth (HPBW)in degrees (7-360 degrees)
theta_3dB_TX = 10; 
% Transmit antenna elevation HPBW in degrees (7-45 degrees)
phi_3dB_TX = 10;
% Receive antenna azimuth HPBW in degrees (7-360 degrees)
theta_3dB_RX = 10; 
% Receive antenna elevation HPBW in degrees (7-45 degrees)
phi_3dB_RX = 10;
% RFBW=0.8;%Bande Passante en GHz
% Tech=1/(RFBW); %en ns
j_complex = sqrt(-1);
RxArrayType = 'ULA'; Nr=1;Wr = Nr; 
dmin = 10; dmax = 100; h_BS=35;
%dans le cas de ULA, ce valeur n'a pas aucun effet sur la matrice du canal
%MIMO mais il a un effet important sur la performance du precodeur
%angulaire
%% simulation parameters
Nr=1;
envType = 'LOS'; 
t = 20;
noise_power=1.3807e-23*RFBW*1e+6*(t+273.15);
save inputParameters 