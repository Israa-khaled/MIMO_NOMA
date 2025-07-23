function getSumRateValues_ADNOMA_2UE_diffUserOrdering_FPA(mu_min,mu_max,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pe,gamma_FPA1)
addpath TempFct
seuil_int=beta_seuil*1000;
for mu=mu_min:mu_max
    [SR_DBS(mu),SR_CB(mu),SR_ZF(mu),SR_2UE_FPA_FCSI(mu),SR_2UE_FPA_ADI(mu),SR_2UE_FPA_dist(mu),SR_2UE_FPA_HW(mu)]= getSRValues_2UE_diffUserOrdering_FPA(mu,Nt,Wt,sceType,TxArrayType,R,f,beta_seuil,scn,Pe,gamma_FPA1)
    save (['struct_ADNOMA_2UE_diffUserOrdering_FPA',sceType,TxArrayType,num2str(Nt),'_',num2str(Wt),'Wt','_',num2str(R),'R',['_',num2str(mu_min),'_',num2str(mu_max),'_',num2str(f),'f',num2str(seuil_int),'_',scn]])
end