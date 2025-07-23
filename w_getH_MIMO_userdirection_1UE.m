
%This function gives us 3 Struct:
%CIR_MIMO_Struct: MIMO matrix for each user 
%AOD_LobePowerSpectrum_Struct: AOD_LobePowerSpectrum for each user
% powerSpectrum(:,1:5)
% AOA_LobePowerSpectrum: is a struct containing Na array (Na =nb of AoA SL), each array contains 
% all the paths which belongs to the same AoA SL.
%AOA_LobePowerSpectrum_Struct: AOA_LobePowerSpectrum for each user 
%powerSpectrum(:,[1:3 6:7]);
%   - powerSpectrum: an array containing the multipath parameters:
%       Column 1: path delays
%       Column 2: path powers, relative to 1 mW
%       Column 3: path phases, in radians
%       Column 4: azimuth angles of departure, in degrees
%       Column 5: zenith angles of departure, in degrees
%       Column 6: azimuth angles of arrival, in degrees
%       Column 7: zenith angles of arrival, in degrees

function[H_ensemble,powerspectrum_struct,TRdist] = w_getH_MIMO_userdirection_1UE(f,RFBW,sceType,envType,dmin,dmax,TXPower,h_BS,Nt,Nr,N,TxArrayType,RxArrayType,Wt,Wr,dTxAnt,dRxAnt,theta_3dB_TX,phi_3dB_TX,theta_3dB_RX,phi_3dB_RX,theta_user,phi_user,TRDistance)

%% Input parameters (subject to change per users' own needs)
freq = num2str(f);
% Barometric Pressure in mbar (1e-5 to 1013.25 mbar)
p = 1013.25; 
% Humidity in % (0-100%)
u = 50; 
% Temperature in degrees Celsius (-100 to 50 degrees Celsius)
t = 20;
% Rain rate in mm/hr (0-150 mm/hr)
RR = 0; 
% Polarization (Co-Pol or X-Pol)
Pol = 'Co-Pol';
% Foliage loss (Yes or No)
Fol = 'No'; 
% Distance within foliage in meters (0-dmin)
dFol = 0; 
% Foliage attenuation in dB/m (0-10 dB/m)
folAtt = 0.4;
%%% Create an output folder 
if exist('NYUSIM_OutputFolder','dir')==0 
    mkdir NYUSIM_OutputFolder
end
%%% Channel Model Parameters
% Free space reference distance in meters
d0 = 1; 
% Speed of light in m/s
c = 3e8;
% Set channel parameters according to the scenario
% UMi LOS
if strcmp(sceType,'UMi') == true && strcmp(envType,'LOS') == true 
% Path loss exponent (PLE)
n = 2; 
% Shadow fading standard deviation in dB
SF = 4.0; 
% Mean angle of departure (AOD)
mu_AOD = 1.9; 
% Mean angle of arrival (AOA)
mu_AOA = 1.8;
% A number between 0 and 1 for generating intra-cluster delays
X_max = 0.2;
% Mean excess delay in ns
mu_tau = 123; 
% Minimum inter-cluster void interval, typically set to 25 ns for outdoor environments
minVoidInterval = 25;
% Per-cluster shadowing in dB
sigmaCluster = 1;
% Time cluster decay constant in ns
Gamma = 25.9; 
% Per-subpath shadowing in dB
sigmaSubpath = 6; 
% Subpath decay constant in ns
gamma = 16.9; 
% Mean zenith angle of departure (ZOD) in degrees
mean_ZOD = -12.6;
% Standard deviation of the ZOD distribution in degrees
sigma_ZOD = 5.9; 
% Standard deviation of the azimuth offset from the lobe centroid
std_AOD_RMSLobeAzimuthSpread = 8.5;
% Standard deviation of the elevation offset from the lobe centroid
std_AOD_RMSLobeElevationSpread = 2.5;
% A string specifying which distribution to use: 'Gaussian' or 'Laplacian'
distributionType_AOD = 'Gaussian'; 
% Mean zenith angle of arrival (ZOA) in degrees
mean_ZOA = 10.8; 
% Standard deviation of the ZOA distribution in degrees
sigma_ZOA = 5.3;
% Standard deviation of the azimuth offset from the lobe centroid
std_AOA_RMSLobeAzimuthSpread = 10.5;
% Standard deviation of the elevation offset from the lobe centroid
std_AOA_RMSLobeElevationSpread = 11.5;
% A string specifying which distribution to use: 'Gaussian' or 'Laplacian'
distributionType_AOA = 'Laplacian';   
% UMi NLOS
elseif strcmp(sceType,'UMi') == true && strcmp(envType,'NLOS') == true
% See the parameter definitions for UMi LOS
n = 3.2; 
SF = 7.0; 
mu_AOD = 1.5; 
mu_AOA = 2.1; 
X_max = 0.5; 
mu_tau = 83;
minVoidInterval = 25; 
sigmaCluster = 3; 
Gamma = 51.0; 
sigmaSubpath = 6;
gamma = 15.5; 
mean_ZOD = -4.9; 
sigma_ZOD = 4.5; 
std_AOD_RMSLobeAzimuthSpread = 11.0;
std_AOD_RMSLobeElevationSpread = 3.0; 
distributionType_AOD = 'Gaussian'; 
mean_ZOA = 3.6; 
sigma_ZOA = 4.8; 
std_AOA_RMSLobeAzimuthSpread = 7.5;
std_AOA_RMSLobeElevationSpread = 6.0; 
distributionType_AOA = 'Laplacian';
% UMa LOS
elseif strcmp(sceType,'UMa') == true && strcmp(envType,'LOS') == true 
% See the parameter definitions for UMi LOS
n = 2; 
SF = 4.0; 
mu_AOD = 1.9; 
mu_AOA = 1.8;
X_max = 0.2; 
mu_tau = 123; 
minVoidInterval = 25; 
sigmaCluster = 1;
Gamma = 25.9; 
sigmaSubpath = 6; 
gamma = 16.9; 
mean_ZOD = -12.6;
sigma_ZOD = 5.9; 
std_AOD_RMSLobeAzimuthSpread = 8.5;
std_AOD_RMSLobeElevationSpread = 2.5;
distributionType_AOD = 'Gaussian'; 
mean_ZOA = 10.8; 
sigma_ZOA = 5.3;
std_AOA_RMSLobeAzimuthSpread = 10.5;
std_AOA_RMSLobeElevationSpread = 11.5;
distributionType_AOA = 'Laplacian'; 
% UMa NLOS
elseif strcmp(sceType,'UMa') == true && strcmp(envType,'NLOS') == true 
% See the parameter definitions for UMi LOS
n = 2.9; 
SF = 7.0; 
mu_AOD = 1.5; 
mu_AOA = 2.1; 
X_max = 0.5; 
mu_tau = 83;
minVoidInterval = 25; 
sigmaCluster = 3; 
Gamma = 51.0; 
sigmaSubpath = 6;
gamma = 15.5; 
mean_ZOD = -4.9; 
sigma_ZOD = 4.5; 
std_AOD_RMSLobeAzimuthSpread = 11.0;
std_AOD_RMSLobeElevationSpread = 3.0; 
distributionType_AOD = 'Gaussian'; 
mean_ZOA = 3.6; 
sigma_ZOA = 4.8; 
std_AOA_RMSLobeAzimuthSpread = 7.5;
std_AOA_RMSLobeElevationSpread = 6.0; 
distributionType_AOA = 'Laplacian';
% RMa LOS
elseif strcmp(sceType,'RMa') == true && strcmp(envType,'LOS') == true
% See the parameter definitions for UMi LOS
SF = 1.7; 
mu_AOD = 1; 
mu_AOA = 1;
X_max = 0.2; 
mu_tau = 123; 
minVoidInterval = 25; 
sigmaCluster = 1;
Gamma = 25.9; 
sigmaSubpath = 6; 
gamma = 16.9; 
mean_ZOD = -12.6;
sigma_ZOD = 5.9; 
std_AOD_RMSLobeAzimuthSpread = 8.5;
std_AOD_RMSLobeElevationSpread = 2.5;
distributionType_AOD = 'Gaussian'; 
mean_ZOA = 10.8; 
sigma_ZOA = 5.3;
std_AOA_RMSLobeAzimuthSpread = 10.5;
std_AOA_RMSLobeElevationSpread = 11.5;
distributionType_AOA = 'Laplacian';
% RMa NLOS
elseif strcmp(sceType,'RMa') == true && strcmp(envType,'NLOS') == true
% See the parameter definitions for UMi LOS
SF = 6.7; 
mu_AOD = 1; 
mu_AOA = 1; 
X_max = 0.5; 
mu_tau = 83;
minVoidInterval = 25; 
sigmaCluster = 3; 
Gamma = 51.0; 
sigmaSubpath = 6;
gamma = 15.5; 
mean_ZOD = -4.9; 
sigma_ZOD = 4.5; 
std_AOD_RMSLobeAzimuthSpread = 11.0;
std_AOD_RMSLobeElevationSpread = 3.0; 
distributionType_AOD = 'Gaussian'; 
mean_ZOA = 3.6; 
sigma_ZOA = 4.8; 
std_AOA_RMSLobeAzimuthSpread = 7.5;
std_AOA_RMSLobeElevationSpread = 6.0; 
distributionType_AOA = 'Laplacian';
end
%% Initialize various settings and parameters
% Structure containing generated CIRs
CIR_SISO_Struct = struct; 
CIR_MIMO_Struct = struct;
% Set plot status
plotStatus = true; 
% Set plot rotation status
plotRotate = false; 
% Determine if spatial plot is needed 
plotSpatial = true;
% Number of multipath components
nPath = zeros(N,1); 
% Best (i.e., smallest) directional path loss
PL_dir_best = zeros(N,1); 
% Directional PDP information
DirPDPInfo = []; 
% Omnidirectional PDP information
OmniPDPInfo = zeros(N,5);
% Run for each RX location, i.e., each channel realization
for CIRIdx = 1:N
clear powerSpectrum PL_dir DirRMSDelaySpread;

    %% Step 1: Generate T-R Separation distance (m) ranging from dmin - dmax.
    %clear TRDistance; TRDistance = getTRSep(dmin,dmax);
    % Set dynamic range, i.e., maximum possible omnidirectional path loss 
    % in dB, according to T-R separation distance. If T-R separation 
    % distance is no larger than 500 m, then set dynamic range as 190 dB, 
    % otherwise set it to 220 dB.
    if TRDistance <= 500
        % Dynamic range in dB
        DR = 190;
    else
        DR = 220;
    end
    % Received power threshod in dBm
    Th = TXPower - DR;
    
    %% Step 2: Generate the total received omnidirectional power (dBm) and 
    % omnidirectional path loss (dB) 
    % non RMa, i.e., UMi or UMa
    if strcmp(sceType,'RMa') == false
    [Pr_dBm, PL_dB]= getRXPower(f,n,SF,TXPower,TRDistance,d0);
    % RMa LOS
    elseif strcmp(sceType,'RMa') == true && strcmp(envType,'LOS') == true 
        PL_dB = 20*log10(4*pi*d0*f*1e9/c) + 23.1*(1-0.03*((h_BS-35)/35))*log10(TRDistance) + SF*randn;
    % RMa NLOS
    elseif strcmp(sceType,'RMa') == true && strcmp(envType,'NLOS') == true 
        PL_dB = 20*log10(4*pi*d0*f*1e9/c) + 30.7*(1-0.049*((h_BS-35)/35))*log10(TRDistance) + SF*randn;
    end
    % Atmospheric attenuation factor
    attenFactor = mpm93_forNYU(f,p,u,t,RR);
    % Path loss incorporating atmospheric attenuation
    PL_dB = getAtmosphericAttenuatedPL(PL_dB,attenFactor,TRDistance);
    % Incorporating cross-polarization
    if strcmp(Pol,'X-Pol') == true
        PL_dB = PL_dB+25;
    end
    % Incorporating foliage loss
    if strcmp(Fol,'Yes') == true
        PL_dB = getFoliageAttenuatedPL(PL_dB,folAtt,dFol);
    end      
    % Calculate received power based on transmit power and path loss
    Pr_dBm = TXPower - PL_dB;
    % Free space path loss
    FSPL = 20*log10(4*pi*d0*f*1e9/c);
    % PLE
    PLE = (PL_dB-FSPL)/(10*log10(TRDistance/d0));
    
    %% Step 3: Generate the number of time clusters N, and number of AOD and
    % AOA spatial lobes
    [numberOfTimeClusters,numberOfAOALobes,numberOfAODLobes] = ...
        getNumClusters_AOA_AOD(mu_AOA,mu_AOD,sceType);
    
    %% Step 4: Generate the number of cluster subpaths M_n for each time cluster
    numberOfClusterSubPaths = getNumberOfClusterSubPaths(numberOfTimeClusters,sceType);
    
    %% Step 5: Generate the intra-cluster subpath delays rho_mn (ns)
    rho_mn = getIntraClusterDelays(numberOfClusterSubPaths,X_max);
    
    %% Step 6: Generate the phases (rad) for each cluster
    phases_mn = getSubpathPhases(rho_mn);
    
    %% Step 7: Generate the cluster excess time delays tau_n (ns)
    tau_n = getClusterExcessTimeDelays(mu_tau,rho_mn,minVoidInterval);
    
    %% Step 8: Generate temporal cluster powers (mW)
    clusterPowers = getClusterPowers(tau_n,Pr_dBm,Gamma,sigmaCluster,Th);
    
    %% Step 9: Generate the cluster subpath powers (mW)
    subpathPowers = ...
        getSubpathPowers(rho_mn,clusterPowers,gamma,sigmaSubpath,envType,Th);
    
    %% Step 10: Recover absolute propagation times t_mn (ns) of each subpath 
    % component
    t_mn = getAbsolutePropTimes(TRDistance,tau_n,rho_mn);
    
    %% Step 11: Recover AODs and AOAs of the multipath components
    [subpath_AODs, cluster_subpath_AODlobe_mapping] = ...
        getSubpathAngles_userdirection_2UEs(numberOfAODLobes,numberOfClusterSubPaths,mean_ZOD,...
        sigma_ZOD,std_AOD_RMSLobeElevationSpread,std_AOD_RMSLobeAzimuthSpread,...
        distributionType_AOD,theta_user);
    [subpath_AOAs, cluster_subpath_AOAlobe_mapping] = ...
        getSubpathAngles(numberOfAOALobes,numberOfClusterSubPaths,mean_ZOA,...
        sigma_ZOA,std_AOA_RMSLobeElevationSpread,std_AOA_RMSLobeAzimuthSpread,...
        distributionType_AOA);
    %% Step 12: Construct the multipath parameters
    powerSpectrumOld = getPowerSpectrum_v1_1(numberOfClusterSubPaths,t_mn,subpathPowers,phases_mn,...
        subpath_AODs,subpath_AOAs,Th);
    %%
    % Adjust power spectrum according to RF bandwidth
    
    %%%%%%%%%% Modification in v 1.6.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % One more output variable "SubpathIndex" is added, which is not being
    % used currently. It is saved for possible future use.
    
    [powerSpectrum,numberOfClusterSubPaths, SubpathIndex] = ...
        getNewPowerSpectrum(powerSpectrumOld,RFBW);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For LOS environment, adjust subpath AoDs and AoAs such that the AoD
    % and AoA of the LOS component are aligned properly
    if strcmp(envType,'LOS') == true
        % Calculate the correct azimuth AoA for LOS component, which should
        % differ from its azimuth AoD by 180 degrees
        % powerSpectrum(1,4) denotes the azimuth AoD for LOS component
        clear correctAzAOA;
        if powerSpectrum(1,4) - 180 > 0
            correctAzAOA = powerSpectrum(1,4) - 180;
        else
            correctAzAOA = powerSpectrum(1,4) + 180;
        end
        % Calculate the difference between the generated azimuth AoA and
        % the correct azimuth AoA
        % powerSpectrum(1,6) is the generated azimuth AoA for LOS component
        clear dAzAOA;
        dAzAOA = powerSpectrum(1,6) - correctAzAOA;
        % Do a global shift of azimuth AoAs
        powerSpectrum(:,6) = powerSpectrum(:,6) - dAzAOA;
        clear azAOA_temp;
        azAOA_temp = powerSpectrum(:,6);
        azAOA_temp(azAOA_temp < 0) = azAOA_temp(azAOA_temp < 0) + 360;
        powerSpectrum(:,6) = azAOA_temp; 
        % Calculate the correct elevation AoA for LOS component, which
        % should be the additive inverse of the corresponding elevation AoD
        clear correctElAOA;
        correctElAOA = -powerSpectrum(1,5); 
        % Calculate the difference between the generated elevation AoA and
        % the correct elevation AoA
        % powerSpectrum(1,7) is the generated elevation AoA for LOS component
        clear dElAOA;
        dElAOA = powerSpectrum(1,7) - correctElAOA;
        % Do a global shift of elevation AoAs
        powerSpectrum(:,7) = powerSpectrum(:,7) - dElAOA;
    end
    
    %% Construct the 3-D lobe power spectra at TX and RX
    %%%%%%%%%% Modification in v 1.6.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In earlier versions, the program crashed when less than 800 MHz
    % bandwidth was specified by the user.
    %
    % v 1.6.1 fixed the problem of previous versions such that the number 
    % of resolvable MPCs decreases as the RF bandwidth is set to be
    % narrower than 800 MHz.
    %
    % However, the v 1.6.1 keeps the spatial resolution (angular information)
    % using 800 MHz RF channel bandwidth despite the user input. 
    % Here PowerSpectrumOld is used for plotting AOAs and AODs of all 
    % resolved MPCs of 800 MHz RF bandwidth regardless of the user bandwidth
    % input in v1.6.1.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AOD_LobePowerSpectrum = getLobePowerSpectrum(numberOfAODLobes,cluster_subpath_AODlobe_mapping,powerSpectrumOld,'AOD');
    AOA_LobePowerSpectrum = getLobePowerSpectrum(numberOfAOALobes,cluster_subpath_AOAlobe_mapping,powerSpectrumOld,'AOA');

    %% Store CIR parameters
    % Multipath delay
    CIR.pathDelays = powerSpectrumOld(:,1);
    % Multipath power
    pathPower = powerSpectrumOld(:,2);
    clear indNaN; indNaN = find(pathPower<=10^(Th/10));
    pathPower(indNaN,:) = 10^(Th/10);
    CIR.pathPowers = pathPower;
    % Multipath phase
    CIR.pathPhases = powerSpectrumOld(:,3);
    % Multipath AOD
    CIR.AODs = powerSpectrumOld(:,4);
    % Multipath ZOD
    CIR.ZODs = powerSpectrumOld(:,5);
    % Multipath AOA
    CIR.AOAs = powerSpectrumOld(:,6);
    % Multipath ZOA
    CIR.ZOAs = powerSpectrumOld(:,7);
    %LOS component
    TRdist(CIRIdx)=TRDistance;
    AOA(CIRIdx)=powerSpectrum(1,6);
    ZOA(CIRIdx)=powerSpectrum(1,7);
    AOD(CIRIdx) = powerSpectrum(1,4);
    ZOD(CIRIdx) = powerSpectrum(1,5);
    % Various global information for this CIR
    % Carrier frequency
    CIR.frequency = freq;
    % Transmit power
    CIR.TXPower = TXPower;
    % Omnidirectional received power in dBm
    CIR.OmniPower = Pr_dBm;
    % Omnidirectional path loss in dB
    CIR.OmniPL = PL_dB;
    % T-R separation distance in meters
    CIR.TRSep = TRDistance;
    % Environment, LOS or NLOS
    CIR.environment = envType;
    % Scenario, UMi, UMa, or RMa
    CIR.scenario = sceType;
    % TX HPBW
    CIR.HPBW_TX = [theta_3dB_TX phi_3dB_TX];
    % RX HPBW
    CIR.HPBW_RX = [theta_3dB_RX phi_3dB_RX];  
    % Store
    CIR_SISO_Struct.(['CIR_SISO_',num2str(CIRIdx)]) = CIR;
    
    %%%%%%%%%%%%%% Modification in v1.6.1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % getLocalCIR_new is modified to generate MIMO channel impulse
    % responses that is stored in CIR_MIMO. H_MIMO stores the CIRs of
    % MPCs resolved at user-specified bandwidth. H is the CIRs of all MPCs
    % resolved at 800 MHz (the widest bandwidth). 
    %
    % In the previous versions of NYUSIM, only the first TX antenna element
    % was used to generate H matrix. In v. 1.6.1, all TX antenna elements 
    % are used, and the H matrix is extended to have a size of Nt x Nr for 
    % each multipath.

    [CIR_MIMO,H,HPowers,HPhases,H_ensemble] = getLocalCIR_new(CIR,powerSpectrumOld,TxArrayType,RxArrayType,Nt,Nr,Wt,Wr,dTxAnt,dRxAnt,RFBW);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CIR_MIMO_Struct.(['CIR_MIMO_',num2str(CIRIdx)]) = CIR_MIMO; H_MIMO = CIR_MIMO.H;
    H_total1(CIRIdx,:)=H_ensemble;
    AOD_LobePowerSpectrum_Struct.(['AOD_LobePowerSpectrum_',num2str(CIRIdx)])=AOD_LobePowerSpectrum;
    AOA_LobePowerSpectrum_Struct.(['AOA_LobePowerSpectrum_',num2str(CIRIdx)])=AOA_LobePowerSpectrum;
    powerspectrum_struct.(['powerspectrum_',num2str(CIRIdx)])=powerSpectrum;
    TR_dist(CIRIdx)=TRDistance;

end 