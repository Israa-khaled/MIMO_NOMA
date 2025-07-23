function getSumRateValues_noma_dbs_OPT1_multiUE_FCSIADI_FPA_FTPA(mu_min,mu_max,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe)
addpath TempFct
seuil_int=beta_seuil*1000;
for mu=mu_min:mu_max
    [sumrate_DBS(mu),sumrate_CB(mu),sumrate_ZF(mu),sumrate_dbsnoma_multiUE_fullCSI(mu),sumrate_dbsnoma_2UE_fullCSI(mu),sumrate_dbsnoma_multiUE_partialCSI(mu),sumrate_dbsnoma_2UE_partialCSI(mu),sumrate_dbsnoma_2UE_FPA1(mu),sumrate_dbsnoma_2UE_FPA2(mu),sumrate_dbsnoma_2UE_FPA3(mu),sumrate_dbsnoma_2UE_FTPA(mu),sumrate_dbsnoma_2UE_FTPA1(mu),sumrate_dbsnoma_2UE_FPA1_ADI(mu),sumrate_dbsnoma_2UE_FPA2_ADI(mu),sumrate_dbsnoma_2UE_FPA3_ADI(mu),sumrate_dbsnoma_2UE_FTPA_ADI(mu),sumrate_dbsnoma_2UE_FTPA_ADI1(mu)]= getsumrate_noma_dbs_OPT1_multiUE_FullPartialCSI_FPA_FTPA(mu,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe);
    save (['struct_nomadbs_OPT1_multiUE_FCSIADI_FPA_FTPA',sceType,TxArrayType,num2str(Nt),'_',num2str(Wt),'Wt','_',num2str(R),'R',['_',num2str(mu_min),'_',num2str(mu_max),'_',num2str(f),'f',num2str(seuil_int),'_',scn]])
end