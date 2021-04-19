function [MU,CPL,RL] = VisCpRL(TTM,PATM,MCLTM)
%This function will firstly calculate density of liquid of single grid.
%Created on 2020-6-17

%TTM: [K] [1,1]
%PATM: absolute pressure [Pa] [1,1]
%MCLTM: temporary MCL of size [1,1]

global MOxide
global Cp
global MV0
global dVdT
global dVdP

%% -------------------- MELT VISCOSITY [Pa.sec] ----------------------

%Normalize oxides to 100 wt%, e.g., 0.541 => 54.1 wt%
All=MCLTM.SiO2+MCLTM.TiO2+MCLTM.Al2O3+MCLTM.FeO+MCLTM.Fe2O3+MCLTM.MnO+MCLTM.MgO+MCLTM.CaO+MCLTM.Na2O+MCLTM.K2O+MCLTM.P2O5+MCLTM.H2O;

SiO2=100.0*MCLTM.SiO2/All;
TiO2=100.0*MCLTM.TiO2/All;
Al2O3=100.0*MCLTM.Al2O3/All;
FeO=100.0*MCLTM.FeO/All;
Fe2O3=100.0*MCLTM.Fe2O3/All;
MnO=100.0*MCLTM.MnO/All;
MgO=100.0*MCLTM.MgO/All;
CaO=100.0*MCLTM.CaO/All;
Na2O=100.0*MCLTM.Na2O/All;
K2O=100.0*MCLTM.K2O/All;
P2O5=100.0*MCLTM.P2O5/All;
H2O=100.0*MCLTM.H2O/All;

%Convert wt% to moles
nSiO2=SiO2/MOxide.SiO2;
nTiO2=TiO2/MOxide.TiO2;
nAl2O3=Al2O3/MOxide.Al2O3;
nFeO=FeO/MOxide.FeO;
nFe2O3=Fe2O3/MOxide.Fe2O3;
nMnO=MnO/MOxide.MnO;
nMgO=MgO/MOxide.MgO;
nCaO=CaO/MOxide.CaO;
nNa2O=Na2O/MOxide.Na2O;
nK2O=K2O/MOxide.K2O;
nP2O5=P2O5/MOxide.P2O5;
nH2O=H2O/MOxide.H2O;

%Mole fraction of species as in following form: [H2O] [K2O, Na2O] [MgO,
%FeO] [CaO, TiO2] [AlO2], i.e., all iron in FeO; all Al in AlO2
nAll=nSiO2+nTiO2+2.0*nAl2O3+nFeO+2.0*nFe2O3+0.0*nMnO+nMgO+nCaO+nNa2O+nK2O+0.0*nP2O5+nH2O;
%all oxides exclude MnO and P2O5
xSiO2=nSiO2/nAll;
xTiO2=nTiO2/nAll;
xAl2O3=2.0*nAl2O3/nAll;
xFeO=nFeO/nAll;
xFe2O3=2.0*nFe2O3/nAll;
xMnO=nMnO/nAll;
xMgO=nMgO/nAll;
xCaO=nCaO/nAll;
xNa2O=nNa2O/nAll;
xK2O=nK2O/nAll;
xP2O5=nP2O5/nAll;
xH2O=nH2O/nAll;

%Slope of each specied group * SiO2 mole fraction * species mole fraction
%Slopes: [AlO2] -> 6.7; [FeO, MgO] -> 3.4; [CaO, TiO2] -> 4.5; [Na2O, K2O]
%-> 2.8; [H2O] -> 2.0

SumSlope=6.7*xSiO2*xAl2O3+...%[AlO2]
    3.4*xSiO2*(xFeO+xFe2O3+xMgO)+...%[FeO, MgO]
    4.5*xSiO2*(xCaO+xTiO2)+...%[CaO, TiO2]
    2.8*xSiO2*(xNa2O+xK2O)+...%[Na2O, K2O]
    2.0*xH2O;%[H2O]
        
MeanSlope=SumSlope/(1.0-xSiO2);

lnLiqVisco=10000.0*MeanSlope/TTM-1.50*MeanSlope-6.40;%T in [K]; LiqVisco in [Poise]
MU=exp(lnLiqVisco)/10.0;%[Pa.sec]   /10.0

%% --------------------- MELT DENSITY [kg/m^3] -------------------------
%Moles of oxides
nRAll=nSiO2+nTiO2+nAl2O3+nFeO+nFe2O3+nMnO+nMgO+nCaO+nNa2O+nK2O+nP2O5+nH2O;%R --> rho
%Molar fraction
xRSiO2=nSiO2/nRAll;
xRTiO2=nTiO2/nRAll;
xRAl2O3=nAl2O3/nRAll;
xRFeO=nFeO/nRAll;
xRFe2O3=nFe2O3/nRAll;
xRMnO=nMnO/nRAll;
xRMgO=nMgO/nRAll;
xRCaO=nCaO/nRAll;
xRNa2O=nNa2O/nRAll;
xRK2O=nK2O/nRAll;
xRP2O5=nP2O5/nRAll;
xRH2O=nH2O/nRAll;

PAG=PATM/10^9;%absolute pressure in GPa

%molar mean volume [m^3/mol]
VL=xRSiO2*(MV0.SiO2+dVdT.SiO2*(TTM-1673.0)+dVdP.SiO2*PAG)+...%SiO2
    xRTiO2*(MV0.TiO2+dVdT.TiO2*(TTM-1673.0)+dVdP.TiO2*PAG)+...%TiO2
    xRAl2O3*(MV0.Al2O3+dVdT.Al2O3*(TTM-1673.0)+dVdP.Al2O3*PAG)+...%Al2O3
    xRFeO*(MV0.FeO+dVdT.FeO*(TTM-1673.0)+dVdP.FeO*PAG)+...%FeO
    xRFe2O3*(MV0.Fe2O3+dVdT.Fe2O3*(TTM-1673.0)+dVdP.Fe2O3*PAG)+...%Fe2O3
    xRMnO*(MV0.MnO+dVdT.MnO*(TTM-1673.0)+dVdP.MnO*PAG)+...%MnO, mostly useless
    xRMgO*(MV0.MgO+dVdT.MgO*(TTM-1673.0)+dVdP.MgO*PAG)+...%MgO
    xRCaO*(MV0.CaO+dVdT.CaO*(TTM-1673.0)+dVdP.CaO*PAG)+...%CaO
    xRNa2O*(MV0.Na2O+dVdT.Na2O*(TTM-1673.0)+dVdP.Na2O*PAG)+...%Na2O
    xRK2O*(MV0.K2O+dVdT.K2O*(TTM-1673.0)+dVdP.K2O*PAG)+...%K2O
    xRP2O5*(MV0.P2O5+dVdT.P2O5*(TTM-1673.0)+dVdP.P2O5*PAG)+...%P2O5, mostly useless
    xRH2O*(MV0.H2O+dVdT.H2O*(TTM-1673.0)+dVdP.H2O*PAG);%H2O
%PAG-1bar ~ PAG

%molar mean mass [g/mol]
ML=xRSiO2*MOxide.SiO2+...
    xRTiO2*MOxide.TiO2+...
    xRAl2O3*MOxide.Al2O3+...
    xRFeO*MOxide.FeO+...
    xRFe2O3*MOxide.Fe2O3+...
    xRMnO*MOxide.MnO+...
    xRMgO*MOxide.MgO+...
    xRCaO*MOxide.CaO+...
    xRNa2O*MOxide.Na2O+...
    xRK2O*MOxide.K2O+...
    xRP2O5*MOxide.P2O5+...
    xRH2O*MOxide.H2O;

%melt density [kg/m^3]
RL=0.001*ML/VL;

%% -------------------- MELT ISOBARIC HEAT CAPACITY [J/kg/K] ----------------------
%Partial Molar Isobaric Heat Capacity for Molten Oxide Components Applicable to
%Silicate Melts at 1 bar. Cp is Approximately Independent of Temperature at T>=1400 K.
%Cp=struct('SiO2',1331.0,'TiO2',1392.0,'Al2O3',1545.0,'Fe2O3',1434.0,'FeO',1100.0,'MnO',1100.0,'MgO',2424.0,'CaO',1781.0,'Na2O',1651.0,'K2O',1030.0,'P2O5',1478.8,'H2O',2278.0);%[J/kg/K]
%NOTE: Data from Table S5.7 in Encyclopedia of Volcanoes, 2015; MnO is
%assumed equal to FeO, and P2O5 is mean of all oxides' Cp except H2O

gfw=(xRSiO2*MOxide.SiO2+xRTiO2*MOxide.TiO2+xRAl2O3*MOxide.Al2O3+xRFe2O3*MOxide.Al2O3+xRFeO*MOxide.FeO+xRMnO*MOxide.MnO+xRMgO*MOxide.MgO+xRCaO*MOxide.CaO+xRNa2O*MOxide.Na2O+xRK2O*MOxide.K2O+xRP2O5*MOxide.P2O5+xRH2O*MOxide.H2O);%[g/mol]

CPL=1000.0*(xRSiO2*Cp.SiO2+xRTiO2*Cp.TiO2+xRAl2O3*Cp.Al2O3+xRFe2O3*Cp.Al2O3+xRFeO*Cp.FeO+xRMnO*Cp.MnO+xRMgO*Cp.MgO+xRCaO*Cp.CaO+xRNa2O*Cp.Na2O+xRK2O*Cp.K2O+xRP2O5*Cp.P2O5+xRH2O*Cp.H2O)/gfw;%[J/kg/K]

%NOTE: An-Di and An-Di-Fo have eutectic temperature > 1270 deg. C (1543 K),
%which means the above CPL is always valid. When An-Di or An-Di-Fo system
%becomes solid, solid CPS have another formula.

end

