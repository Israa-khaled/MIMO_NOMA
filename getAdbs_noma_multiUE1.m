
function [a_BS_noma, azimu]= getAdbs_noma_multiUE1(Clusters_array,azimuthAOD_USER,eleAOD_USER,nb_2UEinCluster,nb_clusters,TxArrayType,Nt,Wt,dTxAnt)
    n_BS = [1:1:Nt]';
    j_complex = sqrt(-1);
    for c=1:nb_2UEinCluster
        azAOD_u=[];phi_u=[];
        for nbUperC=1:length(Clusters_array{c})
            azAOD_u=[azAOD_u,azimuthAOD_USER(Clusters_array{c}(nbUperC))];
            phi_u=[phi_u,eleAOD_USER(Clusters_array{c}(nbUperC))];
        end
        azSort=sort(azAOD_u);
        phiSort=sort(phi_u);
        azAOD_beam=(azSort(1)+azSort(length(Clusters_array{c})))/2;
        azimu(c)=azAOD_beam;
        phi_beam=(phiSort(1)+phiSort(length(Clusters_array{c})))/2;
        if TxArrayType == 'ULA'
            for n_bs = 1:length(n_BS)
                a_BS_noma(n_bs,c) = exp(j_complex.*(n_BS(n_bs)-1).*2.*pi.*dTxAnt.*cos(azAOD_beam*pi/180));    
            end
        elseif TxArrayType == 'URA'
            for n_bs = 1:length(n_BS)
                n_BS_temp=n_BS(n_bs)-1;
                a_BS_noma(n_bs,c) = exp(j_complex.*2.*pi.*dTxAnt.*(abs(mod(n_BS_temp,Wt)).*cos(azAOD_beam*pi/180).*cos(phi_beam*pi/180)+...
                        abs(fix(n_BS_temp/Wt)).*sin(phi_beam*pi/180)));
            end
        end
    end 
    for c=nb_2UEinCluster+1:nb_clusters
        aziAOD=azimuthAOD_USER(Clusters_array{c});
        azimu(c)=aziAOD;
        eleAOD=eleAOD_USER(Clusters_array{c});
        if TxArrayType == 'ULA'
            for n_bs = 1:length(n_BS)
              a_BS_noma(n_bs,c) = exp(j_complex.*(n_BS(n_bs)-1).*2.*pi.*dTxAnt.*cos(aziAOD*pi/180));    
            end
        elseif TxArrayType == 'URA'
            for n_bs = 1:length(n_BS)
                n_BS_temp=n_BS(n_bs)-1;
                a_BS_noma(n_bs,c) = exp(j_complex.*2.*pi.*dTxAnt.*(abs(mod(n_BS_temp,Wt)).*cos(aziAOD*pi/180).*cos(eleAOD*pi/180)+...
                        abs(fix(n_BS_temp/Wt)).*sin(eleAOD*pi/180)));
            end
        end
    end   
end

    
    
    
  