%This script is used to calculate total thermal conductivity (lattice + radiation) of solid
%mineral up to 2000 K, from method of Hofmeister 1999
%Created on 2020-6-21

T=[300:10:2000.0]';%[K]
P=0.0;%[GPa]
a0=0.3407e-4;%thermal expansivity from Fei 1995
a1=0.8674e-8;
a2=-0.7545;

Tintg=a0*(T-298)+0.5*a1*(T.^2-298^2)+a2*(-1.0./T+1.0/298.0);%Eq.(10)
K298=5.2;%7.7%Data from Hofmeister 1999 Table 1
a=0.45;%0.33;
dKTdP=5.2;%4.0;
KT=183;%127.9;%isothermal bulk modulus
ga=1.25;%Gruneisen parameter
KL=K298*(298.0./T).^a.*exp(-Tintg*(1.0/3.0+4.0*ga))*(1.0+dKTdP*P/KT)^((4.0*ga+1.0/3.0)/dKTdP);%lattice only
KR=8.5e-11*T.^3;%radiation only
KLR=KL+KR;%lattice+radiation
plot(T,KLR);