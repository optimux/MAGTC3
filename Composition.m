%This script is used to calculate mass of FeO, Fe2O3 from Table 1 of
%Kuritani 2004
%created 2020-3-31

%SiO2, TiO2, Al2O3, Fe2O3, MnO, MgO, CaO, Na2O, K2O, P2O5, H2O
M=[60.0855 79.867 101.963 159.69 70.938 40.305 56.078 61.9795 94.1966 141.9475 18.016];%Atomic weight [g/mole]

%SiO2, TiO2, Al2O3, Fe2O3, MnO, MgO, CaO, Na2O, K2O, P2O5, H2O
Oxides=[51.4 1.37 18.0 8.59 0.14 5.79 9.72 4.15 0.58 0.28 4.00];%Chemical composition [%]

Mole=Oxides./M;%moles of oxides

SUMM=sum(Mole)+Mole(4);
MF=Mole/(SUMM);%Total iron as FeO
MF(4)=2.0*MF(4);

a=0.218130;
b=13184.7;
c=-4.49933;
dSiO2=-2.15036;
dAl2O3=-8.35163;
dFeOT=-4.49508;
dMgO=-5.43639;
dCaO=0.073113;
dNa2O=3.54148;
dK2O=4.18688;

T=1117.0;%[C]

lnfO2=-2.4648*10^4/(T+273.15)+8.7559;%[Ni-NiO buffer]

lnFe3Fe2=a*lnfO2+b/(T+273.15)+c+dSiO2*MF(1)+dAl2O3*MF(3)+dFeOT*MF(4)+dMgO*MF(6)+dCaO*MF(7)+dNa2O*MF(8)+dK2O*MF(9);%MF(Fe2O3)/MF(FeO)

MFeO=2.0*Mole(4)/(exp(lnFe3Fe2)*2.0+1.0);%mole of FeO
MFe2O3=MFeO*exp(lnFe3Fe2);

FeO=MFeO*71.845;%FeO mass
Fe2O3=MFe2O3*159.69;%Fe2O3 mass