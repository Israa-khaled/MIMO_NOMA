function getSumRateValues_noma_dbs_OPT1_multiUE_FullPartialCSI(mu_min,mu_max,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe)
addpath TempFct
seuil_int=beta_seuil*1000;
for mu=mu_min:mu_max
    [sumrate_DBS(mu),sumrate_CB(mu),sumrate_ZF(mu),sumrate_dbsnoma_multiUE(mu),sumrate_dbsnoma_2UE(mu)]= getsumrate_noma_dbs_OPT1_multiUE_FullPartialCSI(mu,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pmin,Pe);
    save (['struct_nomadbs_OPT1_multiUE_FullPartialCSI',sceType,TxArrayType,num2str(Nt),'_',num2str(Wt),'Wt','_',num2str(R),'R',['_',num2str(mu_min),'_',num2str(mu_max),'_',num2str(f),'f',num2str(seuil_int),'_',scn]], 'sumrate_DBS', 'sumrate_CB','sumrate_ZF','sumrate_dbsnoma_multiUE','sumrate_dbsnoma_2UE','Pe','Pmin')
end