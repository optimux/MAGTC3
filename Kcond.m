%This script is used to calculate thermal conductivity of melt, from Eq.(11) in Ni
%2015
T=[1200:10:1800];%[K]
k=1.380649e-23;%[J/K]
Z=10;
NA=6.022e23;%molecule/mol
V0SiO2=26.86e-6;%m^3/mol
V0MgO=11.69e-6;
V0CaO=16.53e-6;

aTSiO2=0.0/V0SiO2;
aTMgO=3.27e-9/V0MgO;
aTCaO=3.74e-9/V0CaO;
aPSiO2=1.89e-6/V0SiO2;
aPMgO=-0.27e-6/V0MgO;
aPCaO=-0.34e-6/V0CaO;

%Assume 1 mol CaMgSi2O6
V=2.0*(V0SiO2+0.0*(T-1673)-1.89e-6*0.0)+...
    (V0MgO+3.27e-9*(T-1673)+0.27e-6*0.0)+...
    (V0CaO+3.74e-9*(T-1673)+0.34e-6*0.0);%[m^3]
rho=0.001*216.5504./V;%[kg/m^3]
CP=(60.0843*2*1331+40.3044*2424+56.0774*1781)/216.5504;%[J/kg/K]
aT=0.0*2.0+1.0*3.27e-9/11.69e-6+1.0*3.74e-9/16.53e-6;%[1/K]
bT=2.0e-9*1.89e-6/26.86e-6-1.0e-9*0.27e-6/11.69e-6-1.0e-9*0.34e-6/16.53e-6;%[1/Pa]
bS=bT-(T*aT^2)./(rho*CP);
cub=nthroot(Z*NA./V,3);

%dRhodP=-rho*(0.5*dVdPSiO2+0.25*dVdPMgO+0.25*dVdPCaO)./V;
%Kconduc=2.7*k*cub.^2./sqrt(rho.*bS);

Kconduc=2.7*k*cub.^2./sqrt(rho./24.2e9);
plot(T,Kconduc);