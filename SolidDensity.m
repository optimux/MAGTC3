function RSNS = SolidDensity(TTM,PATM)
%This function is used to calculate solid mineral end-members' density.
%Method from Eq.(2-4) in Schutt&Lesher 2006.
%Data from Table 3 in Schutt&Lesher 2006 and Table 4 in
%Korenaga&Korenage2016
%NOTE: all end-members in Schutt&Lesher2006 and Korenaga&Korenaga2016 can
%be calculated for given T and P

%Created on 2020-6-19

%TTM: temperature in K [1,1]
%PATM: absolute pressure in Pa [1,1]

global RS0
global dKdT
global dKdP
global KT0

RSNS=struct('OL',0.0,...%olivine density [kg/m^3]
    'OPX',0.0,...%opx density [kg/m^3]
    'CPX',RS0.Di,...%cpx density [kg/m^3]
    'PL',RS0.An,...%pl density [kg/m^3]
    'ILM',0.0);%ilmentite density [kg/m^3]

TR=298.0;%standard T [K]
PAG=PATM/10^9;%absolute pressure [GPa]
P0G=10^5/10^9;%ambient pressure [GPa]

%NOTE: Schutt&Lesher provides KS0, so convert KS0 to KT0 first at 298 K,
%1 bar: KS0=KT0(1+alpha*gamma*T)
% Fo=128.8/(1+1.15*(0.285e-4+298*1.008e-8-0.384/298^2)*298); --> KT0=127.61
% Fa=138.0/(1+1.12*(0.2386e-4+298*1.153e-8-0.0518/298^2)*298); --> KT0=136.78
% Di=110.5/(1+1.0*(0.232e-4+298*1.88e-8)*298); --> KT0=109.56
% Hd=119.0/(1+1.5*(0.232e-4+298*1.88e-8)*298); --> KT0=117.49

TDi=0.232e-4*(TTM-TR)+0.5*1.88e-8*(TTM^2-TR^2)-0.0*(1.0/TTM-1.0/TR)+0.2*0.0*(TTM^5-TR^5);%Temperature integration of Di
TAn=(2.44e-5-3.1e-7*1.0+1.8e-9*1.0^2)*(TTM-TR)+2.0*(9.0e-9-4.0e-11*1.0)*(0.5*(TTM^2-TR^2)-298.0*(TTM-TR));%Temperature integration of An, from Korenaga2016
PDi=log((KT0.Di+dKdP.Di*(PAG-P0G)+dKdT.Di*(TTM-TR))/(KT0.Di+dKdT.Di*(TTM-TR)));%Pressure integration of Di
PAn=log((KT0.An+dKdP.An*(PAG-P0G)+dKdT.An*(TTM-TR))/(KT0.An+dKdT.An*(TTM-TR)));%Pressure integration of An
RSNS.CPX=RS0.Di*exp(PDi/dKdP.Di-TDi);
RSNS.PL=RS0.An*exp(PAn/dKdP.An-TAn);

RSNS.OL=0.0;
RSNS.OPX=0.0;
RSNS.ILM=0.0;
end

