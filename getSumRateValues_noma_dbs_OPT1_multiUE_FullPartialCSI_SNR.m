function getSumRateValues_noma_dbs_OPT1_multiUE_FullPartialCSI_SNR(SNR_min,SNR_max,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe,MU)
addpath TempFct
seuil_int=beta_seuil*1000;
snr=SNR_min:SNR_max;
for index=1:length(snr)
    index
    a=snr(index)/10;
    NP=Pe/10^a;
    [sumrate_DBS(index),sumrate_CB(index),sumrate_ZF(index),sumrate_dbsnoma_multiUE_fullCSI(index),sumrate_dbsnoma_2UE_fullCSI(index),sumrate_dbsnoma_multiUE_partialCSI(index),sumrate_dbsnoma_2UE_partialCSI(index)]= getsumrate_noma_dbs_OPT1_multiUE_FullPartialCSI_SNR(MU,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe,NP);
    save (['struct_nomadbs_OPT1_multiUE_FullPartialCSI_SNR',sceType,TxArrayType,num2str(Nt),'_',num2str(Wt),'Wt','_',num2str(R),'R',['_',num2str(SNR_min),'_',num2str(SNR_max),'_',num2str(f),'f',num2str(seuil_int),'_',scn]])
end