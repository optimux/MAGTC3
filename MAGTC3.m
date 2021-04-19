%This script is used to calculate An-Di(80) system cooling from above, a
%very similar case in which postulated LMO solidification process may fit
%Created on 2020-6-16
%You're afraid to lose, afraid to try, afraid to change, that's why you're
%going to fail to have great career!

%vx at bottom is suspended, so it's slippery.

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================

clear all

%% ========================= INITIALIZATION ===============================
global NIX
global NIY
global x
global y
global dxI
global dyI
global g
global CE
global DL
global dx
global dy
global RL0
global CL0
global dtb
global FSE
global FSELOG
global P0
global NFSS
global NFSE
global MCL
global MCS
global TCL
global TCS
global TW
global K0
global FRCSMj
global FRCSTe
global Majors
global Minors
global KP
global MOxide
global Cp
global MV0
global dVdT
global dVdP
global RS0
global dKdT
global dKdP
global KT0
global FRCPS
global FSCR1
global FSCR2
global RS
global CPS
global dFSLH
global dFS
global FK
global KS
global FS
global Ea
global Va
global DL0
global HS
global TE
global RSN
global MMins
global MCSN
global MCSEM
global RLE
global CPLE
%global RSECPXPL
global RSE
global CPSE
global FSO
global RSO
global MISS
global TISS
global KST
global TL0

tic;

W=20.0;%width [m]
H=20.0;%height [m]
NIX=40;%number of unit, x-axis
NIY=40;%number of unit, y-axis
dxI=W/NIX;%x-axis interval [m]
dyI=H/NIY;%y-axis interval [m]

%QUICK scheme: grid peclet number Peg=rho*v*dx/mu <=8/3, if rho=3000 kg/m^3,
%v=10^-5 m/sec, mu=0.001 Pa.sec, then dx<=0.09 m, we set dx=0.05 m as maximum
%dx. Tentatively, dx is in 0.01 m order.

dx(1:NIX)=dxI;%[m]
dy(1:NIY)=dyI;%[m]

%Control point (main grid) x-axis location [m]
x=zeros(NIX+2,1);
x(1)=0.0;
x(2)=0.5*dxI;
x(3:NIX+1)=x(2)+[1:NIX-1]*dxI;
x(NIX+2)=x(NIX+1)+0.5*dxI;

%Control point (main grid) y-axis location [m]
y=zeros(NIY+2,1);
y(1)=0.0;
y(2)=0.5*dyI;
y(3:NIY+1)=y(2)+[1:NIY-1]*dyI;
y(NIY+2)=y(NIY+1)+0.5*dyI;

%Initial major elements in liquid [1]
%NOTE: mass fraction (zhi-liang-bai-fen-nong-du), like 51.4 g/100 g=51.4 wt.%=0.514, thus the followings have unit of [1]
%Sample: Di(CaMgSi2O6) 80wt% + An(CaAl2Si2O8) 20wt%
SiO2=0.5303;%SiO2 mass fraction
%NOTE: Pure CaMgSi2O6 has SiO2 0.5549
TiO2=0.0;%TiO2 mass fraction
Al2O3=0.0733;%Al2O3 mass fraction
FeO=0.0;%Fe2O3 mass fraction
Fe2O3=0.0;%Fe2O3 mass fraction
MnO=0.0;%MnO mass fraction
MgO=0.14888;%MgO mass fraction
%NOTE: Pure CaMgSi2O6 has MgO 0.1861
CL0=MgO;
CaO=0.24752;%CaO mass fraction
%NOTE: Pure CaMgSi2O6 has CaO 0.259
Na2O=0.0;%Na2O mass fraction
K2O=0.0;%K2O mass fraction
P2O5=0.0;%P2O5 mass fraction
H2O=0.0;%H2O mass fraction
CE=struct('SiO2',0.5056,'TiO2',0.0,'Al2O3',0.147,'FeO',0.0,'Fe2O3',0.0,'MnO',0.0,'MgO',0.1115,'CaO',0.2359,'Na2O',0.0,'K2O',0.0,'P2O5',0.0,'H2O',0.0);%eutectic composition from fitting of Bowen 1915 (40.1 wt% An) [1]
%NOTE: eutectic temperature is 1271+273.15 K, fitting eutectic MgO concentration in liquid is 11.1479 wt%

%Initial Trace elements in liquid [1]
Sm=1.0e-6;
Nd=1.0e-6;

T0=1357.0+273.15;%Initial T =1630 K > 1623 (TL) [K]
TL0=-1.3351*(CL0*100.0)^2+55.878*CL0*100.0+813.9991+273.15;%liquidus of initial composition [K]
TW=250.0;%lunar cosmological background [K]
TS=1271.0+273.15;%solidus [K]
TE=-1.3351*(100.0*CE.MgO)^2+55.878*CE.MgO*100.0+813.9991+273.15;%eutectic temperature 1544.15 [K]
FSE=zeros(NIY+2,NIX+2);%solid fraction when eutectic
FSELOG=ones(NIY+2,NIX+2);%log of FSE
DS=0.0;%diffusion coefficient in solid matter [m^2/sec]
g=1.62;%lunar gravity [m/sec^2]
PW=1.0e4;%primordial atmosphere pressure [Pa]
K0=5.56e99;%permeability coefficient from Barboza 1998 [m^2]
FSCR1=0.0;%0.34;%total FS at first critical rheology
%VERY IMPORTANT NOTE: If any solid is assumed in place, just set FSCR1<0.0.
%See also NASV.m
FSCR2=0.55;%total FS at second critical rheology

NS=500;%total steps
t=zeros(NS+1,1);%all time sequence [sec]
dt=zeros(NS,1);%all time intervals between successive step [sec]

%Species diffusivity pre-exponential factor from Table s5.10 of Lesher&Spera2015 [m^2/sec]
DL0=struct('SiO2',exp(-12.8),...
    'TiO2',exp(-12.8),...
    'Al2O3',exp(-12.8),...
    'Fe2O3',exp(-4.8),...
    'FeO',exp(-11.2),...
    'MnO',exp(-11.2),...%Lesher2015 + Brady1995 = exp(-3.8)
    'MgO',exp(-9.8),...
    'CaO',exp(-4.272),...%Brady1995, basalt, 1260-1440 celsus degree, 0.1 MPa
    'Na2O',exp(-4.02),...%Brady1995, basalt, 1300-1400 celsus degree, 0.1 MPa
    'K2O',exp(-9.3),...
    'P2O5',exp(-12.652),...%Brady1995, rhyolite, 1200-1500 celsus degree, 800MPa
    'H2O',exp(-9),...
    'Sm',exp(-11.7),...
    'Nd',exp(-13.7));%Lesher2015

%Activation energy of liquid species used for liquid diffusivity from Table s5.10 of Lesher&Spera2015 + Chakraborty2008 [J/mol]
Ea=struct('SiO2',1.67e5,...%Lesher2015
    'TiO2',1.65e5,...
    'Al2O3',3.0e5,...%Van Orman et al., 2001
    'Fe2O3',2.65e5,...%Lesher2015
    'FeO',2.0e5,...%Dohmen&Chakraborty2007
    'MnO',1.66e5,...%Lesher2015 + Brady1995 = 2.01e5
    'MgO',2.0e5,...%Dohmen&Chakraborty2007
    'CaO',1.841e5,...%Brady1995, basalt, 1260-1440 celsus degree, 0.1 MPa
    'Na2O',1.64e5,...%Lesher2015 + Brady1995, basalt, 1300-1400 celsus degree, 0.1 MPa
    'K2O',1.65e5,...
    'P2O5',6.009e5,...%Brady1995, rhyolite, 1200-1500 celsus degree, 800MPa
    'H2O',1.65e5,...
    'Sm',8.4e5,...%Cherniak et al., 1997
    'Nd',1.98e5);%Lesher2015

%Activation volume of liquid species used for liquid diffusivity from Table s5.10 of Lesher&Spera2015 [m^3/mol]
Va=struct('SiO2',-9.0e-6,...
    'TiO2',7.0e-6,...
    'Al2O3',3.05e-6,...%Béjina2003
    'Fe2O3',7.0e-6,...
    'FeO',7.0e-6,...
    'MnO',7.0e-6,...
    'MgO',7.0e-6,...
    'CaO',7.0e-6,...
    'Na2O',7.0e-6,...
    'K2O',7.0e-6,...
    'P2O5',7.0e-6,...
    'H2O',7.0e-6,...
    'Sm',7.0e-6,...
    'Nd',7.0e-6);

%Diffusion coefficient of all species in liquid [m^2/sec]
DL=struct('SiO2',DL0.SiO2*ones(NIY+2,NIX+2),...
    'TiO2',DL0.TiO2*ones(NIY+2,NIX+2),...
    'Al2O3',DL0.Al2O3*ones(NIY+2,NIX+2),...
    'Fe2O3',DL0.Fe2O3*ones(NIY+2,NIX+2),...
    'FeO',DL0.FeO*ones(NIY+2,NIX+2),...
    'MnO',DL0.MnO*ones(NIY+2,NIX+2),...
    'MgO',DL0.MgO*ones(NIY+2,NIX+2),...
    'CaO',DL0.CaO*ones(NIY+2,NIX+2),...
    'Na2O',DL0.Na2O*ones(NIY+2,NIX+2),...
    'K2O',DL0.K2O*ones(NIY+2,NIX+2),...
    'P2O5',DL0.P2O5*ones(NIY+2,NIX+2),...
    'H2O',DL0.H2O*ones(NIY+2,NIX+2),...
    'Sm',DL0.Sm*ones(NIY+2,NIX+2),...
    'Nd',DL0.Nd*ones(NIY+2,NIX+2));

%Molar mass of oxides [g/mol]
MOxide=struct('SiO2',60.0843,...
    'TiO2',79.8658,...
    'Al2O3',101.961276,...
    'Fe2O3',159.6922,...
    'FeO',71.8464,...
    'MnO',70.937449,...
    'MgO',40.3044,...
    'CaO',56.0774,...
    'Na2O',61.97894,...
    'K2O',94.196,...
    'P2O5',141.9446,...
    'H2O',18.01528);

%Molar mass of mineral end-members [g/mol]
MMins=struct('Fo',140.6931,...%Forsterite   Mg2SiO4
    'Fa',203.7771,...%Fayalite    Fe2SiO4
    'Di',216.5504,...%Diopside   CaMgSi2O6
    'Hd',248.0924,...%Hedenbergite   CaFeSi2O6
    'En',200.7774,...%Enstatite   Mg2Si2O6
    'Fs',263.8614,...%Ferrosilite    Fe2Si2O6
    'Py',403.1274,...%Pyrope    Mg3Al2[SiO4]3
    'Gs',450.4464,...%Grossular   Ca3Al2[SiO4]3
    'Sp',142.2657,...%Spinel   MgAl2O4
    'Us',223.5586,...%Ulvospinel   TiFe(II)2O4
    'An',278.2073,...%Anorthite   CaAl2Si2O8
    'Ab',262.2230);%Albite   NaAlSi3O8
%NOTE: more end-members are in alphaMELTS forum/mannual, Schutt&Lesher 2006 and Korenaga 2016

%Partial Molar Isobaric Heat Capacity for Molten Oxide Components Applicable to
%Silicate Melts at 1 bar. Cp is Approximately Independent of Temperature at T>=1400 K [J/kg/K]
%Cp in [J/kg/K]
% Cp=struct('SiO2',1331.0,...
%     'TiO2',1392.0,...
%     'Al2O3',1545.0,...
%     'Fe2O3',1434.0,...
%     'FeO',1100.0,...
%     'MnO',1100.0,...
%     'MgO',2424.0,...
%     'CaO',1781.0,...
%     'Na2O',1651.0,...
%     'K2O',1030.0,...
%     'P2O5',1478.8,...
%     'H2O',2278.0);

%Cp in [J/gfw/K]
Cp=struct('SiO2',80.0,...%[J/gfw/K]
    'TiO2',111.8,...
    'Al2O3',157.6,...
    'Fe2O3',229.0,...
    'FeO',78.9,...
    'MnO',78.9,...
    'MgO',99.7,...
    'CaO',99.9,...
    'Na2O',102.3,...
    'K2O',97.0,...
    'P2O5',113.51,...
    'H2O',41.0);
%NOTE: Data from Table S5.7 in Encyclopedia of Volcanoes, 2015; MnO is
%assumed equal to FeO, and P2O5 is mean of all oxides' Cp except H2O

%Partial molar volumes of melten oxides at 1673 K [m^3/mol]
MV0=struct('SiO2',26.86e-6,...
    'TiO2',23.16e-6,...
    'Al2O3',37.42e-6,...
    'Fe2O3',42.13e-6,...
    'FeO',13.65e-6,...
    'MnO',13.65e-6,...
    'MgO',11.69e-6,...
    'CaO',16.53e-6,...
    'Na2O',28.88e-6,...
    'K2O',45.07e-6,...
    'P2O5',25.904e-6,...
    'H2O',26.27e-6);
%Isobaric expansivity of melten oxides dV/dT [m^3/mol/K]
dVdT=struct('SiO2',0.0e-9,...
    'TiO2',7.24e-9,...
    'Al2O3',0.0e-9,...
    'Fe2O3',9.09e-9,...
    'FeO',2.92e-9,...
    'MnO',2.92e-9,...
    'MgO',3.27e-9,...
    'CaO',3.74e-9,...
    'Na2O',7.68e-9,...
    'K2O',12.08e-9,...
    'P2O5',4.89e-9,...
    'H2O',9.46e-9);
%Isothermal compressibility of melten oxides dV/dP [m^3/mol/GPa]
dVdP=struct('SiO2',-1.89e-6,...
    'TiO2',-2.31e-6,...
    'Al2O3',-2.26e-6,...
    'Fe2O3',-2.53e-6,...
    'FeO',-0.45e-6,...
    'MnO',-0.45e-6,...
    'MgO',0.27e-6,...
    'CaO',0.34e-6,...
    'Na2O',-2.4e-6,...
    'K2O',-6.75e-6,...
    'P2O5',-1.84e-6,...
    'H2O',-3.15e-6);
%NOTE: Data from Table S5.1 in Encyclopedia of Volcanoes, 2015; MnO is
%assumed equal to FeO, and P2O5 is mean of all oxides' volumes except H2O

%Solid mineral end-members density at 298 K, 1 bar. Data from Schutt&Lesher2006 and Korenaga2016 [kg/m^3]
RS0=struct('An',2765.0,'Ab',2615.0,'Fo',3230.5,'Fa',4391.5,'Di',3277.0,'Hd',3656.0);

%Solid mineral end-members Isothermal bulk modulus at 298 K, 1 bar. Data from
%Schutt&Lesher2006 and Korenaga2016 [GPa]; some data is also available in
%Linda 2008
%NOTE: Schutt&Lesher provides KS0, so convert KS0 to KT0 first at 298 K,
%1 bar: KS0=KT0(1+alpha*gamma*T)
% Fo=128.8/(1+1.15*(0.285e-4+298*1.008e-8-0.384/298^2)*298); --> 127.61
% Fa=138.0/(1+1.12*(0.2386e-4+298*1.153e-8-0.0518/298^2)*298); --> 136.78
% Di=110.5/(1+1.0*(0.232e-4+298*1.88e-8)*298); --> 109.56
% Hd=119.0/(1+1.5*(0.232e-4+298*1.88e-8)*298); --> 117.49
KT0=struct('An',81.6,'Ab',57.6,'Fo',127.61,'Fa',136.78,'Di',109.56,'Hd',117.49);

%Solid mineral end-members isothermal bulk modulus pressure derivative dKT/dP [GPa/GPa]
dKdP=struct('An',4.0,'Ab',4.0,'Fo',4.2,'Fa',4.7,'Di',4.8,'Hd',4.0);
%NOTE: second derivative with respect to pressure is neglected.

%Solid mineral end-members isothermal bulk modulus temperature derivative dKT/dT [GPa/K]
dKdT=struct('An',-0.02,'Ab',-0.02,'Fo',-0.017,'Fa',-0.0204,'Di',-0.0205,'Hd',-0.0205);

Majors={'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO','CaO','Na2O','K2O','P2O5','H2O'};%Majors oxides names
Minors={'Sm','Nd'};%trace elements
Mins={'OL','OPX','CPX','PL','ILM'};%minerals

ChemMajor=zeros(length(Majors),1);%mass fraction for major oxides
ChemMinor=zeros(length(Minors),1);%mass fraction for trace elements

alloxides=0.0;
for i=1:length(Majors)
    ChemMajor(i)=eval(Majors{i});
    alloxides=alloxides+ChemMajor(i);
end

if(abs(alloxides-1.0)>=0.01)
    fprintf('The sum of all oxides is not 100!\n');
end

for i=1:length(Minors)
    ChemMinor(i)=eval(Minors{i});
end

%======================= Main field variables =============================
%T=T(t,x,y);
%CL=CL(t,x,y);
%CS=CS(t,x,y);%CS not time dependent in solid
%FL=FL(t,x,y);
%FS=FS(t,x,y);
%dFS=dFS(t,x,y)
%VX=VX(t,x,y);
%VY=VY(t,x,y);
%P=P(t,x,y);

%RD==RecorD

%------------------------- Temperature Field ------------------------------
T=T0*ones(NIY+2,NIX+2);%T in chamber [K]
%NOTE 1: if top boundary is in physical contact with cold rocks, and heat is transfered
%by heat conduction k*dT/dz
%T(1,2:NIX+1)=2.0*TW-T(2,2:NIX+1);%numerical top cool wall [K]
%NOTE 2: if top boundary is in contact with cold atmosphere, then heat is carried away
%by air convection, h*dT where h is convective heat transfer coefficient, h is 5~25 W/m^2/K
%T(1,2:NIX+1)=TW;%top cool wall [K]

T(1,2:NIX+1)=T(2,2:NIX+1);%psuedo top cool wall [K]
T(NIY+2,2:NIX+1)=T(NIY+1,2:NIX+1);%bottom insulated wall [K]
T(2:NIY+1,1)=T(2:NIY+1,2);%left insulate [K]
T(2:NIY+1,NIX+2)=T(2:NIY+1,NIX+1);%right insulate [K]
T(1,1)=T(1,2);
T(1,NIX+2)=T(1,NIX+1);
T(NIY+2,1)=T(NIY+2,2);
T(NIY+2,NIX+2)=T(NIY+2,NIX+1);
TRD=zeros(NIY+2,NIX+2,NS+1);%T records of all steps [K]
TRD(:,:,1)=T;%1st T record [K]

%----------------------- Major Elements in Liquid -------------------------
MCL=struct('SiO2',SiO2*ones(NIY+2,NIX+2),...%SiO2 in liquid [1]
    'TiO2',TiO2*ones(NIY+2,NIX+2),...%TiO2 in liquid [1]
    'Al2O3',Al2O3*ones(NIY+2,NIX+2),...%Al2O3 in liquid [1]
    'FeO',FeO*ones(NIY+2,NIX+2),...%FeO in liquid [1]
    'Fe2O3',Fe2O3*ones(NIY+2,NIX+2),...%Fe2O3 in liquid [1]
    'MnO',MnO*ones(NIY+2,NIX+2),...%MnO in liquid [1]
    'MgO',MgO*ones(NIY+2,NIX+2),...%MgO in liquid [1]
    'CaO',CaO*ones(NIY+2,NIX+2),...%CaO in liquid [1]
    'Na2O',Na2O*ones(NIY+2,NIX+2),...%Na2O in liquid [1]
    'K2O',K2O*ones(NIY+2,NIX+2),...%K2O in liquid [1]
    'P2O5',P2O5*ones(NIY+2,NIX+2),...%P2O5 in liquid [1]
    'H2O',H2O*ones(NIY+2,NIX+2));%H2O in liquid [1]

MCLRD=struct('SiO2',SiO2*ones(NIY+2,NIX+2,NS+1),...%SiO2 record in liquid [1]
    'TiO2',TiO2*ones(NIY+2,NIX+2,NS+1),...%TiO2 record in liquid [1]
    'Al2O3',Al2O3*ones(NIY+2,NIX+2,NS+1),...%Al2O3 record in liquid [1]
    'FeO',FeO*ones(NIY+2,NIX+2,NS+1),...%FeO record in liquid [1]
    'Fe2O3',Fe2O3*ones(NIY+2,NIX+2,NS+1),...%Fe2O3 record in liquid [1]
    'MnO',MnO*ones(NIY+2,NIX+2,NS+1),...%MnO record in liquid [1]
    'MgO',MgO*ones(NIY+2,NIX+2,NS+1),...%MgO record in liquid [1]
    'CaO',CaO*ones(NIY+2,NIX+2,NS+1),...%CaO record in liquid [1]
    'Na2O',Na2O*ones(NIY+2,NIX+2,NS+1),...%Na2O record in liquid [1]
    'K2O',K2O*ones(NIY+2,NIX+2,NS+1),...%K2O record in liquid [1]
    'P2O5',P2O5*ones(NIY+2,NIX+2,NS+1),...%P2O5 record in liquid [1]
    'H2O',H2O*ones(NIY+2,NIX+2,NS+1));%H2O record in liquid [1]

%----------------------- Trace Elements in Liquid -------------------------
TCL=struct('Sm',Sm*ones(NIY+2,NIX+2),...%Sm in liquid [1]
    'Nd',Nd*ones(NIY+2,NIX+2));...%Nd in liquid [1]
    
TCLRD=struct('Sm',Sm*ones(NIY+2,NIX+2,NS+1),...%Sm record in liquid [1]
    'Nd',Nd*ones(NIY+2,NIX+2,NS+1));...%Nd record in liquid [1]
    
%----------------------- Major Elements in Solid --------------------------
%It should be noted that MCS.OXIDE.ANY is based on mass weight average.
MCS=struct('SiO2',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%SiO2 in OL, OPX, CPX, PL, ILM [1]
    'TiO2',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%TiO2 in OL, OPX, CPX, PL, ILM [1]
    'Al2O3',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%Al2O3 in OL, OPX, CPX, PL, ILM [1]
    'FeO',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%FeO in OL, OPX, CPX, PL, ILM [1]
    'Fe2O3',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%Fe2O3 in OL, OPX, CPX, PL, ILM [1]
    'MnO',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%MnO in OL, OPX, CPX, PL, ILM [1]
    'MgO',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%MgO in OL, OPX, CPX, PL, ILM [1]
    'CaO',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%CaO in OL, OPX, CPX, PL, ILM [1]
    'Na2O',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%Na2O in OL, OPX, CPX, PL, ILM [1]
    'K2O',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%K2O in OL, OPX, CPX, PL, ILM [1]
    'P2O5',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%P2O5 in OL, OPX, CPX, PL, ILM [1]
    'H2O',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)));%H2O in OL, OPX, CPX, PL, ILM [1]

MCSN=MCS;%New MCS

%MCS in End-Members initialization
MCSEM=struct('SiO2',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'TiO2',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'Al2O3',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'FeO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'Fe2O3',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'MnO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'MgO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'CaO',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'Na2O',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'K2O',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'P2O5',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0),...
    'H2O',struct('Fo',0.0,'Fa',0.0,'Di',0.0,'Hd',0.0,'En',0.0,'Fs',0.0,'Py',0.0,'Gs',0.0,'Sp',0.0,'Us',0.0,'An',0.0,'Ab',0.0));
%Calculate MCS in End-Members
OxideEndMember;

MCSRD=struct('SiO2',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%SiO2 record in OL, OPX, CPX, PL, ILM [1]
    'TiO2',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%TiO2 record in OL, OPX, CPX, PL, ILM [1]
    'Al2O3',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%Al2O3 record in OL, OPX, CPX, PL, ILM [1]
    'FeO',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%FeO record in OL, OPX, CPX, PL, ILM [1]
    'Fe2O3',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%Fe2O3 record in OL, OPX, CPX, PL, ILM [1]
    'MnO',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%MnO record in OL, OPX, CPX, PL, ILM [1]
    'MgO',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%MgO record in OL, OPX, CPX, PL, ILM [1]
    'CaO',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%CaO record in OL, OPX, CPX, PL, ILM [1]
    'Na2O',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%Na2O record in OL, OPX, CPX, PL, ILM [1]
    'K2O',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%K2O record in OL, OPX, CPX, PL, ILM [1]
    'P2O5',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%P2O5 record in OL, OPX, CPX, PL, ILM [1]
    'H2O',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)));%H2O record in OL, OPX, CPX, PL, ILM [1]

%----------------------- Trace Elements in Solid --------------------------
TCS=struct('Sm',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)),...%Sm in OL, OPX, CPX, PL, ILM [1]
    'Nd',struct('OL',zeros(NIY+2,NIX+2),'OPX',zeros(NIY+2,NIX+2),'CPX',zeros(NIY+2,NIX+2),'PL',zeros(NIY+2,NIX+2),'ILM',zeros(NIY+2,NIX+2)));...%Nd in OL, OPX, CPX, PL, ILM [1]
    
TCSRD=struct('Sm',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)),...%Sm record in OL, OPX, CPX, PL, ILM [1]
    'Nd',struct('OL',zeros(NIY+2,NIX+2,NS+1),'OPX',zeros(NIY+2,NIX+2,NS+1),'CPX',zeros(NIY+2,NIX+2,NS+1),'PL',zeros(NIY+2,NIX+2,NS+1),'ILM',zeros(NIY+2,NIX+2,NS+1)));%Nd record in OL, OPX, CPX, PL, ILM [1]

%Major Residual Internal Storage Species including solid and liquid [kg/m^3]
MISS=struct('SiO2',zeros(NIY+2,NIX+2),...
    'TiO2',zeros(NIY+2,NIX+2),...
    'Al2O3',zeros(NIY+2,NIX+2),...
    'FeO',zeros(NIY+2,NIX+2),...
    'Fe2O3',zeros(NIY+2,NIX+2),...
    'MnO',zeros(NIY+2,NIX+2),...
    'MgO',zeros(NIY+2,NIX+2),...
    'CaO',zeros(NIY+2,NIX+2),...
    'Na2O',zeros(NIY+2,NIX+2),...
    'K2O',zeros(NIY+2,NIX+2),...
    'P2O5',zeros(NIY+2,NIX+2),...
    'H2O',zeros(NIY+2,NIX+2));
%Initialization after FL, RL, CL.

%Trace Residual Internal Storage Species including solid and liquid [kg/m^3]
TISS=struct('Sm',zeros(NIY+2,NIX+2),...
    'Nd',zeros(NIY+2,NIX+2));
%Initialization after FL, RL, CL.

%------------------------ Volume Fraction ---------------------------------
FL=ones(NIY+2,NIX+2);%liquid volume fraction [1]
FL(1,:)=0.0;%real top boundary
FLRD=zeros(NIY+2,NIX+2,NS+1);%FL records of all steps [1]
FLRD(:,:,1)=FL;%1st record [1]

FS=struct('OL',zeros(NIY+2,NIX+2),...%olivine volume fraction at current step [1]
    'OPX',zeros(NIY+2,NIX+2),...%opx volume fraction at current step [1]
    'CPX',zeros(NIY+2,NIX+2),...%cpx volume fraction at current step [1]
    'PL',zeros(NIY+2,NIX+2),...%pl volume fraction at current step [1]
    'ILM',zeros(NIY+2,NIX+2));%ilmentite volume fraction at current step [1]

FS.OL(1,:)=0.0;%real top boundary
FS.OPX(1,:)=0.0;
FS.CPX(1,:)=(80.0/RS0.Di)/(80.0/RS0.Di+20.0/RS0.An);%Di 80 wt% system, RS0.is at 298 K, 1 bar
FS.PL(1,:)=(20.0/RS0.An)/(80.0/RS0.Di+20.0/RS0.An);%Di 80 wt% system, RS0.is at 298 K, 1 bar
FS.ILM(1,:)=0.0;

FSRD=struct('OL',zeros(NIY+2,NIX+2,NS+1),...%olivine volume fraction record [1]
    'OPX',zeros(NIY+2,NIX+2,NS+1),...%opx volume fraction record [1]
    'CPX',zeros(NIY+2,NIX+2,NS+1),...%cpx volume fraction record [1]
    'PL',zeros(NIY+2,NIX+2,NS+1),...%pl volume fraction record [1]
    'ILM',zeros(NIY+2,NIX+2,NS+1));%ilmentite volume fraction record [1]

%Old FS for FindNSE, to find when minerals begin to form and end
FSO=struct('OL',zeros(NIY+2,NIX+2),...%olivine volume fraction at old step [1]
    'OPX',zeros(NIY+2,NIX+2),...%opx volume fraction at old step [1]
    'CPX',zeros(NIY+2,NIX+2),...%cpx volume fraction at old step [1]
    'PL',zeros(NIY+2,NIX+2),...%pl volume fraction at old step [1]
    'ILM',zeros(NIY+2,NIX+2));%ilmentite volume fraction at old step [1]

dFSLH=struct('OL',zeros(NIY+2,NIX+2),...%olivine volume fraction increment related to latent heat [1]
    'OPX',zeros(NIY+2,NIX+2),...%opx volume fraction increment related to latent heat [1]
    'CPX',zeros(NIY+2,NIX+2),...%cpx volume fraction increment related to latent heat [1]
    'PL',zeros(NIY+2,NIX+2),...%pl volume fraction increment related to latent heat [1]
    'ILM',zeros(NIY+2,NIX+2));%ilmentite volume fraction increment related to latent heat [1]

dFS=struct('OL',zeros(NIY+2,NIX+2),...%total (=dFSSM+dFSLH) olivine volume fraction increment [1]
    'OPX',zeros(NIY+2,NIX+2),...%total (=dFSSM+dFSLH) opx volume fraction increment [1]
    'CPX',zeros(NIY+2,NIX+2),...%total (=dFSSM+dFSLH) cpx volume fraction increment [1]
    'PL',zeros(NIY+2,NIX+2),...%total (=dFSSM+dFSLH) pl volume fraction increment [1]
    'ILM',zeros(NIY+2,NIX+2));%total (=dFSSM+dFSLH) ilmentite volume fraction increment [1]
dFSRD=struct('OL',zeros(NIY+2,NIX+2,NS+1),...%total olivine volume fraction increment [1]
    'OPX',zeros(NIY+2,NIX+2,NS+1),...%total opx volume fraction increment [1]
    'CPX',zeros(NIY+2,NIX+2,NS+1),...%total cpx volume fraction increment [1]
    'PL',zeros(NIY+2,NIX+2,NS+1),...%total pl volume fraction increment [1]
    'ILM',zeros(NIY+2,NIX+2,NS+1));%total ilmentite volume fraction increment [1]

%-------------------------- Momentum Field --------------------------------
PA=PW*ones(NIY+2,NIX+2);%pressure in liquid [Pa]; A==Absolute
%NOTE: Pressure initialization after density initialization

CPL=zeros(NIY+2,NIX+2);%liquid specific heat capacity [J/kg/K]
RL=zeros(NIY+2,NIX+2);%density of liquid [kg/m^3]
MUO=zeros(NIY+2,NIX+2);%Old Dynamic viscosity of liquid [Pa.sec]
%Calculate dynamic viscosity MUO [Pa.sec], liquid specific heat capacity CPL [J/kg/K] and
%liquid density RL [kg/m^3] at given T [K], PA [Pa] and MCL [1]
[MUO,CPL,RL]=VisCpRLs(T,PA,MCL,1.0-FL);
%IMPORTANT NOTE: Temperature should be changed according to temperature boundary condition; velocity should be changed according to FREE, NO SLIP or impermeable
%boundary condition; RL, RS, MU, FS, FL, dFS, PD, KS, KL, CPL, CPS, MCL, MCS, TCL, TCS, RP in ghost cells are set equal to inner layer.

MURD=zeros(NIY+2,NIX+2,NS+1);
MURD(:,:,1)=MUO;

%---------- Pressure Initialization ----------
%NOTE: It should be noted that for a deep magma, say 100s km, pressure effect on denisty can not be neglected!!!
PA(2,:)=0.5*RL(2,:).*FL(2,:)*g*dy(1)+PW;%first row in chamber pressure [Pa]; A==absolute pressure
PA(1,:)=2.0*PW-PA(2,:);%this determines liquid surface is exactly PW Pa [Pa]
for i=1:NIX+2
    for j=3:NIY+1
        PA(j,i)=PA(j-1,i)+0.5*g*dy(j-2)*RL(j-1,i)*FL(j-1,i)+0.5*g*dy(j-1)*RL(j,i)*FL(j,i);%Assuming no solid at t=0
    end
end
PA(NIY+2,:)=PA(NIY+1,:)+0.5*g*dy(NIY)*RL(NIY+1,:)+0.5*g*dy(NIY)*RL(NIY+2,:);%set dy(NIY+1)==dy(NIY)
%NOTE: this makes P at real physical boundary equal to
%P(NIY+1,:)+0.5*g*dy(NIY)*RL(NIY+1,:) since cells (NIY+1,:) are identical
%to cells (NIY+2,:)
P0=PA;
PD=PA-PA;%Dynamics (or relative) pressure in liquid with respect to initial pressure [Pa]
PDRD=zeros(NIY+2,NIX+2,NS+1);%Dynamic (or relative) P records of all steps [Pa]
PDRD(:,:,1)=PD;%1st record [Pa]

%MUOE=MUO;%Effective dynamic viscosity of mixture [Pa.sec]

[~,CPLE,RLE]=VisCpRL(TE,PW,CE);
%eutectic CPL [J/kg/K]
% CPLE=(0.401*2.0*MOxide.SiO2/MMins.An+0.599*2.0*MOxide.SiO2/MMins.Di)*Cp.SiO2+...
%     (0.401*MOxide.CaO/MMins.An+0.599*MOxide.CaO/MMins.Di)*Cp.CaO+...
%     (0.599*MOxide.MgO/MMins.Di+0.0)*Cp.MgO+...
%     (0.0+0.401*MOxide.Al2O3/MMins.An)*Cp.Al2O3;

VX=zeros(NIY+2,NIX+1);%liquid x-axis velocity [m/sec]
VX(1,:)=-VX(2,:);%Top NO SLIP boundary
VX(NIY+2,:)=VX(NIY+1,:);%Bottom FREE boundary
%VERY IMPORTANT NOTE: the surface of LMO, if no crystals yet, has a FREE boundary condition for VX, while bottom part of LMO has a NO SLIP boundary condition
% VX(1,:)=VX(2,:);%Top FREE boundary
% VX(NIY+2,:)=VX(NIY+1,:);%Bottom FREE boundary
VXRD=zeros(NIY+2,NIX+1,NS+1);%VX records of all steps [m/sec]
VXRD(:,:,1)=VX;%1st record [m/sec]

VY=zeros(NIY+1,NIX+2);%liquid y-axis velocity [m/sec]
VY(:,1)=VY(:,2);%left FREE boundary
VY(:,NIX+2)=VY(:,NIX+1);%right FREE boundary
% VY(:,1)=-VY(:,2);%left NO SLIP boundary
% VY(:,NIX+2)=-VY(:,NIX+1);%right NO SLIP boundary
VYRD=zeros(NIY+1,NIX+2,NS+1);%VX records of all steps [m/sec]
VYRD(:,:,1)=VY;%1st record [m/sec]

%---------------------- Minerals mass left in liquid --------------------------
ML=struct('OL',zeros(NIY+2,NIX+2),...%OL mass left at current step [kg]
    'OPX',zeros(NIY+2,NIX+2),...%opx mass left at current step [kg]
    'CPX',0.8*RL.*ones(NIY+2,NIX+2),...%cpx mass left at current step [kg]
    'PL',0.2*RL.*ones(NIY+2,NIX+2),...%pl mass left at current step [kg]
    'ILM',zeros(NIY+2,NIX+2));%ilm mass left at current step [kg]

MLF=struct('OL',zeros(NIY+2,NIX+2),...%OL mass fraction left at current step [kg]
    'OPX',zeros(NIY+2,NIX+2),...%opx mass fraction left at current step [kg]
    'CPX',0.8*ones(NIY+2,NIX+2),...%cpx mass fraction left at current step [kg]
    'PL',0.2*ones(NIY+2,NIX+2),...%pl mass fraction left at current step [kg]
    'ILM',zeros(NIY+2,NIX+2));%ilm mass fraction left at current step [kg]
%======================= Main field variables =============================

%======================= Auxillary field variables ========================

%--------------------------- Load GeoInfo ---------------------------------
% %perl('run_alphamelts.command','-f','F:\HT\Kuritani\KuritaniSetting.txt','-p','F:\HT\Kuritani\','-b','F:\HT\Kuritani\KuritaniBatch.txt');
% %GeoInfo;%write alphaMELTS data as MAT file
% %load GeoConst;%load necessary database, like partition coefficient
%
% %Load Composition Evolution info from alphaMELTS
% melt=load('melt.mat');%liquid: Temperature mass viscosity SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O Mg#
% OL=load('OL.mat');%olivine: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O
% OPX=load('OPX.mat');%orthopyroxene: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O
% CPX=load('CPX.mat');%clinopyroxene: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O
% PL=load('PL.mat');%feldspar: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O
% ILM=load('ILM.mat');%spinel: Temperature mass SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 H2O

%Trace element solid/liquid partition coefficient, from Elkins-Tanton 2008
KPSmOL=0.0006;%OL
KPSmOPX=0.02;%OPX
KPSmCPX=0.3;%CPX
KPSmPL=0.03;%PL
KPSmILM=0.001;%ILM

KPNdOL=0.0006;%OL
KPNdOPX=0.009;%OPX
KPNdCPX=0.2;%CPX
KPNdPL=0.04;%PL
KPNdILM=0.001;%ILM

KP=struct('Nd',struct('OL',KPNdOL*ones(NIY+2,NIX+2),'OPX',KPNdOPX*ones(NIY+2,NIX+2),'CPX',KPNdCPX*ones(NIY+2,NIX+2),'PL',KPNdPL*ones(NIY+2,NIX+2),'ILM',KPNdILM*ones(NIY+2,NIX+2)),...%Nd partition coefficient of OL, OPX, CPX, PL, ILM
    'Sm',struct('OL',KPSmOL*ones(NIY+2,NIX+2),'OPX',KPSmOPX*ones(NIY+2,NIX+2),'CPX',KPSmCPX*ones(NIY+2,NIX+2),'PL',KPSmPL*ones(NIY+2,NIX+2),'ILM',KPSmILM*ones(NIY+2,NIX+2)));%Sm partition coefficient of OL, OPX, CPX, PL, ILM

%CS are valid only if FS>0, otherwise CS should be NaN, however NaN can not
%be used in the following routine, so we use CL*KP
for i=1:NIX+2
    for j=1:NIY+2
        TCS.Sm.OL(j,i)=TCL.Sm(j,i)*KP.Sm.OL(j,i);
        TCS.Sm.OPX(j,i)=TCL.Sm(j,i)*KP.Sm.OPX(j,i);
        TCS.Sm.CPX(j,i)=TCL.Sm(j,i)*KP.Sm.CPX(j,i);
        TCS.Sm.PL(j,i)=TCL.Sm(j,i)*KP.Sm.PL(j,i);
        TCS.Sm.ILM(j,i)=TCL.Sm(j,i)*KP.Sm.ILM(j,i);
        
        TCS.Nd.OL(j,i)=TCL.Nd(j,i)*KP.Nd.OL(j,i);
        TCS.Nd.OPX(j,i)=TCL.Nd(j,i)*KP.Nd.OPX(j,i);
        TCS.Nd.CPX(j,i)=TCL.Nd(j,i)*KP.Nd.CPX(j,i);
        TCS.Nd.PL(j,i)=TCL.Nd(j,i)*KP.Nd.PL(j,i);
        TCS.Nd.ILM(j,i)=TCL.Nd(j,i)*KP.Nd.ILM(j,i);
    end
end

%------------------------- Medium properties ------------------------------
%KS(T,CS)=KS(t,x,y)
%KL(T,CL)=KL(t,x,y)
%CPS(T,CS)=CPS(t,x,y)
%CPL(T,CL)=CPL(t,x,y)
%RS(T,CS)=RS(t,x,y)
%RL(T,CL)=RL(t,x,y)
%RM=RM(t,x,y)
%CM=CM(t,x,y)
%LH(CS)=LH(t,x,y)
%KP(CL)=KP(t,x,y)
%TL(CL)=TL(t,x,y)
%KFL(FL)=KFL(t,x,y)

%NOTE: these properties generally are functions of temperature; when ghost
%cell is used, temperature in ghost cells may lead to different properties
%than main domain.

CPS=struct('OL',zeros(NIY+2,NIX+2),...%olivine specific heat capacity [J/kg/K]
    'OPX',zeros(NIY+2,NIX+2),...%opx specific heat capacity [J/kg/K]
    'CPX',zeros(NIY+2,NIX+2),...%cpx specific heat capacity [J/kg/K]
    'PL',zeros(NIY+2,NIX+2),...%pl specific heat capacity [J/kg/K]
    'ILM',zeros(NIY+2,NIX+2));%ilmentite specific heat capacity [J/kg/K]

CPSE=struct('OL',0.0,...%eutectic olivine specific heat capacity [J/kg/K]
    'OPX',0.0,...%eutectic opx specific heat capacity [J/kg/K]
    'CPX',0.0,...%eutectic cpx specific heat capacity [J/kg/K]
    'PL',0.0,...%eutectic pl specific heat capacity [J/kg/K]
    'ILM',0.0);%eutectic ilmentite specific heat capacity [J/kg/K]

%NOTE 1: CpS Data from Constitution of the Moon: 5. Constraints on composition,
%density, temperature, and radius of a core, by Kuskov 1998
%NOTE 2: T in Kelvin
%NOTE 3: CPS at ghost cells are set equal to inner layer since the
%medium surrouding this cavity is not given, hence CPS is not given by
%T in ghost cells, we use Border.m.
CPS.PL=290.90+0.0276*T-34080000.0./T.^2+5218000000.0./T.^3+29625.0./T;%[J/mol/K]
CPS.PL(1,:)=290.90+0.0276*TW-34080000.0/TW^2+5218000000.0/TW^3+29625.0/TW;%top boundary [J/mol/K]
CPS.PL=CPS.PL/0.2782073;%CaAl2Si2O8: 278.2073 g/mol [J/kg/K]
CPS.CPX=157.3+0.0000205*T-1372950.0./T.^2-1010.5./sqrt(T);%[J/mol/K]
CPS.CPX(1,:)=157.3+0.0000205*TW-1372950.0/TW^2-1010.5/sqrt(TW);%[J/mol/K]
CPS.CPX=CPS.CPX/0.1082752;%Ca0.5Mg0.5SiO3: 108.2752 g/mol [J/kg/K]

CPSE.PL=290.90+0.0276*TE-34080000.0/TE^2+5218000000.0/TE^3+29625.0/TE;%[J/mol/K]
CPSE.PL=CPSE.PL/0.2782073;%CaAl2Si2O8: 278.2073 g/mol [J/kg/K]
CPSE.CPX=157.3+0.0000205*TE-1372950.0/TE^2-1010.5/sqrt(TE);%[J/mol/K]
CPSE.CPX=CPSE.CPX/0.1082752;%Ca0.5Mg0.5SiO3: 108.2752 g/mol [J/kg/K]

RS=SolidDensitys(T,PA);%solid phase density [kg/m^3]
% RST=SolidDensity(TW,PW);%TW boundary solid density [kg/m^3]
% RS.OL(1,:)=RST.OL;%top boundary
% RS.OPX(1,:)=RST.OPX;
% RS.CPX(1,:)=RST.CPX;
% RS.PL(1,:)=RST.PL;
% RS.ILM(1,:)=RST.ILM;
%NOTE 1: RS is a struct
% RS=struct('OL',zeros(NIY+2,NIX+2),...%olivine density [kg/m^3]
%     'OPX',zeros(NIY+2,NIX+2),...%opx density [kg/m^3]
%     'CPX',RS0.Di*ones(NIY+2,NIX+2),...%cpx density [kg/m^3]
%     'PL',RS0.An*ones(NIY+2,NIX+2),...%pl density [kg/m^3]
%     'ILM',zeros(NIY+2,NIX+2));%ilmentite density [kg/m^3]
%NOTE 2: RS at ghost cells are set equal to inner layer since the
%medium surrouding this cavity is not given, hence RS is not given by
%T in ghost cells, we use Border.m in SolidDensitys.m

RSN=RS;
RSE=SolidDensity(TE,PW);
%RSE=struct('CPX',0.0,'PL',0.0);

%RSE=1.0/(0.599/RSE.CPX+0.401/RSE.PL);%eutectic mean solid density [kg/m^3]

RSRD=struct('OL',zeros(NIY+2,NIX+2,NS+1),...%olivine density [kg/m^3]
    'OPX',zeros(NIY+2,NIX+2,NS+1),...%opx density [kg/m^3]
    'CPX',RS0.Di*ones(NIY+2,NIX+2,NS+1),...%cpx density [kg/m^3]
    'PL',RS0.An*ones(NIY+2,NIX+2,NS+1),...%pl density [kg/m^3]
    'ILM',zeros(NIY+2,NIX+2,NS+1));%ilmentite density [kg/m^3]

RSRD.OL(:,:,1)=RS.OL;%1st record [kg/m^3]
RSRD.OPX(:,:,1)=RS.OPX;%1st record [kg/m^3]
RSRD.CPX(:,:,1)=RS.CPX;%1st record [kg/m^3]
RSRD.PL(:,:,1)=RS.PL;%1st record [kg/m^3]
RSRD.ILM(:,:,1)=RS.ILM;%1st record [kg/m^3]
RSB=RS;%solid density for Boussinesq approximation

FR=zeros(NIY+2,NIX+2);%FS*RS --> total solid density at each cell [kg/m^3]
FR=RS.OL.*FS.OL+RS.OPX.*FS.OPX+RS.CPX.*FS.CPX+RS.PL.*FS.PL+RS.ILM.*FS.ILM;
FRRD=zeros(NIY+2,NIX+2,NS+1);%record total solid density at each cell [kg/m^3]
FRRD(:,:,1)=FR;

FRO=FR;

RLRD=zeros(NIY+2,NIX+2,NS+1);%liquid density records of all steps [kg/m^3]
RLRD(:,:,1)=RL;%1st record [kg/m^3]
RLB=RL;%liquid density for Boussinesq approximation

RL0=RL;%liquid reference density at t=0 sec

RM=RL.*FL+RS.OL.*FS.OL+RS.OPX.*FS.OPX+RS.CPX.*FS.CPX+RS.PL.*FS.PL+RS.ILM.*FS.ILM;%solid + liquid mean density of each cell [kg/m^3]
RMRD=zeros(NIY+2,NIX+2,NS+1);%mean density records of all steps [kg/m^3]
RMRD(:,:,1)=RM;%1st record

KS=Conducts(T,PA);
KST=Conduct(TW,PW);
    KS.OL(1,:)=KST.OL;
    KS.OPX(1,:)=KST.OPX;
    KS.CPX(1,:)=KST.CPX;
    KS.PL(1,:)=KST.PL;
    KS.ILM(1,:)=KST.ILM;
%NOTE 1: KS is a struct
% KS=struct('OL',zeros(NIY+2,NIX+2),...%olivine thermal conductivity [W/m/K]
%     'OPX',zeros(NIY+2,NIX+2),...%opx thermal conductivity [W/m/K]
%     'CPX',1.4*ones(NIY+2,NIX+2),...%cpx thermal conductivity [W/m/K]
%     'PL',1.4*ones(NIY+2,NIX+2),...%pl thermal conductivity [W/m/K]
%     'ILM',zeros(NIY+2,NIX+2));%ilmentite thermal conductivity [W/m/K]
%NOTE 2: KS at ghost cells are set equal to inner layer since the
%medium surrouding this cavity is not given, hence KS is not given by
%T in ghost cells, we use Border.m in SolidConduct.m.
%NOTE 3: surprisingly, ilmenite has a very low thermal conductivity,
%see 'Thermal conductivity of titanium slags' by Heimo 2019
%The following data is from Heimo 2019 T1 sample:
%Temperature [celsius] Total thermal conductivity [W/m/K]
%26.93602694	2.195054945
%101.010101	  2.087912088
%200.3367003	1.989010989
%299.6632997	1.865384615
%398.989899	  1.758241758
%500	1.684065934
%599.3265993	1.700549451
%700.3367003	1.758241758
%799.6632997	1.832417582
%900.6734007	1.857142857
%1000	1.881868132
%1101.010101	1.947802198
%The above data may be fitted by 3rd polynomial, like:
%k=-2.728e-12*T^3 +5.345499e-9*T^2-1.1817e-6*T+2.225736; T in celsius, k in W/m/K

FK=struct('OL',zeros(NIY+2,NIX+2),...%olivine FS*KS [W/m/K]
    'OPX',zeros(NIY+2,NIX+2),...%opx FS*KS [W/m/K]
    'CPX',zeros(NIY+2,NIX+2),...%cpx FS*KS [W/m/K]
    'PL',zeros(NIY+2,NIX+2),...%pl FS*KS [W/m/K]
    'ILM',zeros(NIY+2,NIX+2));%ilmentite FS*KS [W/m/K]

FK.OL=FS.OL.*KS.OL;
FK.OPX=FS.OPX.*KS.OPX;
FK.CPX=FS.CPX.*KS.CPX;
FK.PL=FS.PL.*KS.PL;
FK.ILM=FS.ILM.*KS.ILM;

FKT=zeros(NIY+2,NIX+2);%FS*KS --> total solid thermal conductivity at each cell [W/m/K]
FKT=FK.OL+FK.OPX+FK.CPX+FK.PL+FK.ILM;
%NOTE: real top boundary has invariant composition and temperature, its thermal conductivity stays constant

mf0=80.0/216.5504/(80.0/216.5504+20.0/278.2073);%Di initial mole fraction
mf=mf0*ones(NIY+2,NIX+2);%Di mole fraction
DAn=zeros(NIY+2,NIX+2);%An thermal diffusivity [m^2/sec]
KL=zeros(NIY+2,NIX+2);%liquid thermal conductivity [W/m/K]
DAn=(0.36+0.4*exp(-(T-273.0)/300.0))*1.0e-6;%[m^2/sec]
KL=(DAn*(1.0-mf)+0.57e-6*mf).*RL.*CPL+8.5e-11*T.^3;%molar fraction averaged thermal conductivity of liquid Klat+Krad[W/m/K]
KL(1,1:NIX+2)=KL(2,1:NIX+2);%numerical top boundary
%KL=(0.8*0.290e-6+0.2*0.358e-6)*RL.*CPL;%mass fraction averaged thermal conductivity of liquid [W/m/K]
%NOTE 1:
%An=0.358e-6 m^2/sec glass remelted lattice thermal diffusivity from Hofmeister 2009 (not valid)
%Di=0.290e-6 m^2/sec glass remelted lattice thermal diffusivity from Hofmeister 2009 (not valid)
%NOTE 2:
%An = Table 2 of Branlund 2012 = Eq.(7) in Christopher J. Grose = (0.36+0.4*exp(-(T-273.0)/300.0))*1.0e-6
%Di=0.57e-6 m^2/sec crystal melted to liquid lattice thermal diffusivity from Hofmeister 2008
%NOTE 3: RL and CPL have been set equal at ghost cells to inner layer by
%VisCpRLs.m, hence KL at ghost cells are equal to inner layer.
%NOTE 4: if thermal conductivity has a theory of k=sum(ki*Mi) where Mi is oxide molar fraction,
%then all thermal properties would be generalizedbased on oxide molar fraction. However, no such
%theory could be applicable, so we have to turn to thermal condutivities of mineral endmembers.
%Since melt can not be deveided mathematically into a series of mineral endmembers , we'd better
%treat it as an entity (RL*CPL). Diffusivity of this entity could be a molar fraction mean diffusivity
%of mineral endmembers.

HS=struct('OL',1.01e6*ones(NIY+2,NIX+2),...%Fo latent heat Lesher&Spera2015 [J/kg]
    'OPX',7.29e5*ones(NIY+2,NIX+2),...%En latent heat Lesher&Spera2015 [J/kg]
    'CPX',6.36e5*ones(NIY+2,NIX+2),...%Di latent heat in Lesher&Spera2015 [J/kg]
    'PL',4.78e5*ones(NIY+2,NIX+2),...%An latent heat Lesher&Spera2015 [J/kg]
    'ILM',1.43e5*ones(NIY+2,NIX+2));%ilm latent heat Lesher&Spera2015 [J/kg]
%NOTE: latent heat is assumed to be pressure independent. There's only one
%example of pressure-dependent latent heat (Indium), see Höhne 1996. On the
%other hand, latent heat itself is not temperature dependent.

%HSE=0.599*HS.CPX(2,2)+0.401*HS.PL(2,2);%eutectic solid heat fusion [J/kg]

KFL=zeros(NIY+2,NIX+2);%permeability [m^2]
% Blake-Kozeny-Carman
for i=1:NIX+2
    for j=1:NIY+2
        Ga=1.0;%(0.5+atan(100.0*(1.0-FSCR1-FLNTM(j,i)))/pi)^-5;%rheology transition factor is 5, Ga should be omitted here since no solid motion.
        KFL(j,i)=Ga*K0*(FL(j,i)^3/(1.0-FL(j,i))^2);
    end
end

%---------------------------- Other variables -----------------------------
%NFSS=step at which solidification begins
%NFSE=step at which solidification ends

%To record step at which solid begins to form
NFSS=struct('OL',ones(NIY+2,NIX+2,'int16'),...
    'OPX',ones(NIY+2,NIX+2,'int16'),...
    'CPX',ones(NIY+2,NIX+2,'int16'),...
    'PL',ones(NIY+2,NIX+2,'int16'),...
    'ILM',ones(NIY+2,NIX+2,'int16'));

%To record step at which no more solid to form
NFSE=struct('OL',ones(NIY+2,NIX+2,'int16'),...
    'OPX',ones(NIY+2,NIX+2,'int16'),...
    'CPX',ones(NIY+2,NIX+2,'int16'),...
    'PL',ones(NIY+2,NIX+2,'int16'),...
    'ILM',ones(NIY+2,NIX+2,'int16'));
%NOTE: these two variables above are used to cope with inhomogeneous
%species distribution in solid

%integration FS*RS*CpS of solid phases [J/m^3/K]
FRCPS=struct('OL',zeros(NIY+2,NIX+2),...%olivine
    'OPX',zeros(NIY+2,NIX+2),...%opx
    'CPX',zeros(NIY+2,NIX+2),...%cpx
    'PL',zeros(NIY+2,NIX+2),...%pl
    'ILM',zeros(NIY+2,NIX+2));%ilmentite
FRCPS.OL=FS.OL.*RS.OL.*CPS.OL;
FRCPS.OPX=FS.OPX.*RS.OPX.*CPS.OPX;
FRCPS.CPX=FS.CPX.*RS.CPX.*CPS.CPX;
FRCPS.PL=FS.PL.*RS.PL.*CPS.PL;
FRCPS.ILM=FS.ILM.*RS.ILM.*CPS.ILM;

%FS*RS*CS of major elements record in solid (SM+LH) at each step
%NOTE: major and trace elements in solid can be tracked; you may add more here
%Mj==Major
FRCSMj=struct('SiO2',zeros(NIY+2,NIX+2),...%SiO2 in total solid [1]
    'TiO2',zeros(NIY+2,NIX+2),...%TiO2 in total solid [1]
    'Al2O3',zeros(NIY+2,NIX+2),...%Al2O3 in total solid [1]
    'FeO',zeros(NIY+2,NIX+2),...%FeO in total solid [1]
    'Fe2O3',zeros(NIY+2,NIX+2),...%Fe2O3 in total solid [1]
    'MnO',zeros(NIY+2,NIX+2),...%MnO in total solid [1]
    'MgO',zeros(NIY+2,NIX+2),...%MgO in total solid [1]
    'CaO',zeros(NIY+2,NIX+2),...%CaO in total solid [1]
    'Na2O',zeros(NIY+2,NIX+2),...%Na2O in total solid [1]
    'K2O',zeros(NIY+2,NIX+2),...%K2O in total solid [1]
    'P2O5',zeros(NIY+2,NIX+2),...%P2O5 in total solid [1]
    'H2O',zeros(NIY+2,NIX+2));...%H2O in total solid [1]
    
%FS*RS*CS of trace elements record in solid (SM+LH) at each step
%NOTE: major and trace elements in solid can be tracked; you may add more here
%Te==Trace elements
FRCSTe=struct('Sm',zeros(NIY+2,NIX+2),...%Sm in total solid [1]
    'Nd',zeros(NIY+2,NIX+2));...%Nd in total solid [1]
    
%Mean solid concentration of major and trace elements in solid at 1->n steps [wt %]
CSM=struct('SiO2',zeros(NIY+2,NIX+2),...%SiO2 in 1->n mean solid [wt%]
    'TiO2',zeros(NIY+2,NIX+2),...%TiO2 in 1->n mean solid [wt%]
    'Al2O3',zeros(NIY+2,NIX+2),...%Al2O3 in 1->n mean solid [wt%]
    'FeO',zeros(NIY+2,NIX+2),...%FeO in 1->n mean solid [wt%]
    'Fe2O3',zeros(NIY+2,NIX+2),...%Fe2O3 in 1->n mean solid [wt%]
    'MnO',zeros(NIY+2,NIX+2),...%MnO in 1->n mean solid [wt%]
    'MgO',zeros(NIY+2,NIX+2),...%MgO in 1->n mean solid [wt%]
    'CaO',zeros(NIY+2,NIX+2),...%CaO in 1->n mean solid [wt%]
    'Na2O',zeros(NIY+2,NIX+2),...%Na2O in 1->n mean solid [wt%]
    'K2O',zeros(NIY+2,NIX+2),...%K2O in 1->n mean solid [wt%]
    'P2O5',zeros(NIY+2,NIX+2),...%P2O5 in 1->n mean solid [wt%]
    'H2O',zeros(NIY+2,NIX+2));...%H2O in 1->n mean solid [wt%]
    
% RSCPS=zeros(NIY+2,NIX+2);%RS*CPS, current step
% RLCPL=zeros(NIY+2,NIX+2);%RL*CPL, current step
RESE=zeros(NIY+2,NIX+2);%residual of energy check
RESM=zeros(NIY+2,NIX+2);%residual of mass check

for i=1:NIX+2
    for j=1:NIY+2
        
        %-------------------------- FORTRAN -------------------------------
        %initial total FS*RS*CS for major elements in solid, but now MCS.ANY is 0.
        FRCSMj.SiO2(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.SiO2.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.SiO2.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.SiO2.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.SiO2.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.SiO2.ILM(j,i);
        FRCSMj.TiO2(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.TiO2.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.TiO2.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.TiO2.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.TiO2.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.TiO2.ILM(j,i);
        FRCSMj.Al2O3(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.Al2O3.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.Al2O3.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.Al2O3.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.Al2O3.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.Al2O3.ILM(j,i);
        FRCSMj.FeO(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.FeO.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.FeO.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.FeO.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.FeO.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.FeO.ILM(j,i);
        FRCSMj.Fe2O3(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.Fe2O3.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.Fe2O3.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.Fe2O3.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.Fe2O3.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.Fe2O3.ILM(j,i);
        FRCSMj.MnO(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.MnO.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.MnO.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.MnO.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.MnO.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.MnO.ILM(j,i);
        FRCSMj.MgO(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.MgO.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.MgO.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.MgO.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.MgO.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.MgO.ILM(j,i);
        FRCSMj.CaO(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.CaO.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.CaO.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.CaO.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.CaO.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.CaO.ILM(j,i);
        FRCSMj.Na2O(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.Na2O.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.Na2O.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.Na2O.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.Na2O.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.Na2O.ILM(j,i);
        FRCSMj.K2O(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.K2O.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.K2O.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.K2O.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.K2O.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.K2O.ILM(j,i);
        FRCSMj.P2O5(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.P2O5.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.P2O5.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.P2O5.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.P2O5.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.P2O5.ILM(j,i);
        FRCSMj.H2O(j,i)=FS.OL(j,i)*RS.OL(j,i)*MCS.H2O.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*MCS.H2O.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*MCS.H2O.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*MCS.H2O.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*MCS.H2O.ILM(j,i);
        
        %initial total FS*RS*CS for minor elements in solid
        FRCSTe.Sm(j,i)=FS.OL(j,i)*RS.OL(j,i)*TCS.Sm.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*TCS.Sm.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*TCS.Sm.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*TCS.Sm.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*TCS.Sm.ILM(j,i);
        FRCSTe.Nd(j,i)=FS.OL(j,i)*RS.OL(j,i)*TCS.Nd.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)*TCS.Nd.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)*TCS.Nd.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)*TCS.Nd.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i)*TCS.Nd.ILM(j,i);
        
        %-------------------------- MATLAB --------------------------------
        %         %FRCS for major elements
        %         for k=1:length(Majors)
        %             Chems=0.0;
        %             for m=1:length(Mins)
        %                 cmd=['FS.',Mins{m},'(j,i)*RS.',Mins{m},'(j,i)*MCS.',Majors{k},'.',Mins{m},'(j,i)'];
        %                 Chems=Chems+eval(cmd);
        %             end
        %             cmd=['FRCSMj.',Majors{k},'(j,i)=Chems'];
        %             eval(cmd);
        %         end
        %
        %         %FS*RS*TCS of trace elements in solid
        %         for k=1:length(Minors)
        %             Chems=0.0;
        %             for m=1:length(Mins)
        %                 cmd=['FS.',Mins{m},'(j,i)*RS.',Mins{m},'(j,i)*TCS.',Minors{k},'.',Mins{m},'(j,i)'];
        %                 Chems=Chems+eval(cmd);
        %             end
        %             cmd=['FRCSTe.',Minors{k},'(j,i)=Chems'];
        %             eval(cmd);
        %         end
        
        %Major Residual Internal Storage Speices [kg/m^3]
        MISS.SiO2(j,i)=FL(j,i)*RL(j,i)*MCL.SiO2(j,i)+FRCSMj.SiO2(j,i);
        MISS.TiO2(j,i)=FL(j,i)*RL(j,i)*MCL.TiO2(j,i)+FRCSMj.TiO2(j,i);
        MISS.Al2O3(j,i)=FL(j,i)*RL(j,i)*MCL.Al2O3(j,i)+FRCSMj.Al2O3(j,i);
        MISS.FeO(j,i)=FL(j,i)*RL(j,i)*MCL.FeO(j,i)+FRCSMj.FeO(j,i);
        MISS.Fe2O3(j,i)=FL(j,i)*RL(j,i)*MCL.Fe2O3(j,i)+FRCSMj.Fe2O3(j,i);
        MISS.MnO(j,i)=FL(j,i)*RL(j,i)*MCL.MnO(j,i)+FRCSMj.MnO(j,i);
        MISS.MgO(j,i)=FL(j,i)*RL(j,i)*MCL.MgO(j,i)+FRCSMj.MgO(j,i);
        MISS.CaO(j,i)=FL(j,i)*RL(j,i)*MCL.CaO(j,i)+FRCSMj.CaO(j,i);
        MISS.Na2O(j,i)=FL(j,i)*RL(j,i)*MCL.Na2O(j,i)+FRCSMj.Na2O(j,i);
        MISS.K2O(j,i)=FL(j,i)*RL(j,i)*MCL.K2O(j,i)+FRCSMj.K2O(j,i);
        MISS.P2O5(j,i)=FL(j,i)*RL(j,i)*MCL.P2O5(j,i)+FRCSMj.P2O5(j,i);
        MISS.H2O(j,i)=FL(j,i)*RL(j,i)*MCL.H2O(j,i)+FRCSMj.H2O(j,i);
        
        %Trace Residual Internal Storage Species [kg/m^3]
        TISS.Sm(j,i)=FL(j,i)*RL(j,i)*TCL.Sm(j,i)+FRCSTe.Sm(j,i);
        TISS.Nd(j,i)=FL(j,i)*RL(j,i)*TCL.Nd(j,i)+FRCSTe.Nd(j,i);
        
    end
end

QMRD=zeros(NS+1,1);
%======================= Auxillary field variables ========================

fprintf('Trace Elements: ');
for i=1:length(Minors)
    fprintf('%s ',Minors{i});
end
fprintf('\n');

% aviobj=VideoWriter('MAGTC33.avi');
% open(aviobj);

%============================ Load Metadata ===============================
%NOTE: load metadata from .mat files when either initiallization or
%debugging

T(1,1:NIX+2)=2.0*TW-T(2,1:NIX+2);%numerical boundary
n0=0;
 
% load('RL.mat');%liquid density, non-structure
% load('RS.mat');%solid density, structure
% load('T0.mat');
% load('FL.mat');%liquid volume fraction, non-structure
% load('FS.mat');%sold volume fraction, structure
% load('dFS.mat');%total soild volume fraction increment, structure
% load('dFSLH.mat');%latent soild volume fraction increment, structure
% load('VX.mat');%non-structure
% load('VY.mat');%non-structure
% load('PD.mat');%dynamic (relative) pressure, non-structure
% load('FRCPS.mat');%structure
% load('FRCSMj.mat');%structure
% load('FRCSTe.mat');%structure
% load('FK.mat');%structure
% load('FR.mat');%non-structure
% load('CPL.mat');%non-structure
% load('CPS.mat');%structure
% load('MCS.mat');%structure
% load('TCS.mat');%structure
% load('KS.mat');%structure
% load('KL.mat');%non-structure
% load('MCL.mat');%non-structure
% load('TCL.mat');%non-structure
% load('MU.mat');%non-structure
% load('MUO.mat');%non-structure
% load('FLO.mat');%non-structure
% load('RLO.mat');%non-structure
% load('FRO.mat');%non-structure
% load('MISS.mat');%structure
% load('TISS.mat');%structure
% load('NFSS.mat');%structure
% load('NFSE.mat');%structure
% load('CSM.mat');%structure
% load('time.mat');%non-structure
% load('DL.mat');%non-structure
% load('RSB.mat');%non-structure
% load('RLB.mat');%non-structure
% n0=length(t);%how many dtb in previous run
% t=[t;zeros(NS,1)];

% load('E:\demo13\20214142026\RL.mat');%liquid density, non-structure
% load('E:\demo13\20214142026\RS.mat');%solid density, structure
% load('E:\demo13\20214142026\T0.mat');
% load('E:\demo13\20214142026\FL.mat');%liquid volume fraction, non-structure
% load('E:\demo13\20214142026\FS.mat');%sold volume fraction, structure
% load('E:\demo13\20214142026\dFS.mat');%total soild volume fraction increment, structure
% load('E:\demo13\20214142026\dFSLH.mat');%latent soild volume fraction increment, structure
% load('E:\demo13\20214142026\VX.mat');%non-structure
% load('E:\demo13\20214142026\VY.mat');%non-structure
% load('E:\demo13\20214142026\PD.mat');%dynamic (relative) pressure, non-structure
% load('E:\demo13\20214142026\FRCPS.mat');%structure
% load('E:\demo13\20214142026\FRCSMj.mat');%structure
% load('E:\demo13\20214142026\FRCSTe.mat');%structure
% load('E:\demo13\20214142026\FK.mat');%structure
% load('E:\demo13\20214142026\FR.mat');%non-structure
% load('E:\demo13\20214142026\CPL.mat');%non-structure
% load('E:\demo13\20214142026\CPS.mat');%structure
% load('E:\demo13\20214142026\MCS.mat');%structure
% load('E:\demo13\20214142026\TCS.mat');%structure
% load('E:\demo13\20214142026\KS.mat');%structure
% load('E:\demo13\20214142026\KL.mat');%non-structure
% load('E:\demo13\20214142026\MCL.mat');%non-structure
% load('E:\demo13\20214142026\TCL.mat');%non-structure
% load('E:\demo13\20214142026\MU.mat');%non-structure
% load('E:\demo13\20214142026\MUO.mat');%non-structure
% load('E:\demo13\20214142026\FLO.mat');%non-structure
% load('E:\demo13\20214142026\RLO.mat');%non-structure
% load('E:\demo13\20214142026\FRO.mat');%non-structure
% load('E:\demo13\20214142026\MISS.mat');%structure
% load('E:\demo13\20214142026\TISS.mat');%structure
% load('E:\demo13\20214142026\NFSS.mat');%structure
% load('E:\demo13\20214142026\NFSE.mat');%structure
% load('E:\demo13\20214142026\CSM.mat');%structure
% load('E:\demo13\20214142026\time.mat');%non-structure
% load('E:\demo13\20214142026\DL.mat');%non-structure
% load('E:\demo13\20214142026\RSB.mat');%non-structure
% load('E:\demo13\20214142026\RLB.mat');%non-structure

%OriginPlot(VX,VY,T);

% delete MCLRD.mat
% delete MCSRD.mat
% delete PDRD.mat
% delete TRD.mat
% delete VYRD.mat
% delete VXRD.mat
% delete FSRD.mat
% delete FLRD.mat
% delete dFSRD.mat
% delete dtb.mat


%% ============================ Main loop =================================
for n=1:NS
    fprintf('Calculating Step: %5d   -->  %5d\n',n-1,n);
    
    %--------------------------- T-FS-CL ----------------------------------
    %Integrate solid properties: intg(FS*RS*CPS), intg(FS*RS*MCS), intg(FS*RS*TCS), intg(FS*KS), intg(FS*RS) [NIY+2,NIX+2]
    [FRCPST,FKT,FR] = INTSP(FR,CPL,RL);
    
    %Estimate time step [sec]
    dtb=TSTEP(FRCPST,FL,FS,RL,RS,CPL,KL,FKT,MUO,VX,VY,n);%VSX,VSY,
    
    %Old properties used in P-V coupling with Non_Boussinesq approximation
    FLO=FL;
    RLO=RL;
    RSO=RS;
    %Old FS for FindNSE to determine when minerals begin to form and end
    FSO=FS;
    
    %Solid Cumulation from Solid Movement
    %dFSSMT=SolidAssemble;%(FL);%total solid cumulation due to movement [NIY,NIX]
    
    dt(n)=dtb;
    t(n0+n+1)=t(n0+n)+dt(n);
    
    %Old mean density [kg/m^3]
    RM=RL.*FL+FR;
    RMRD(:,:,n+1)=RM;
    
    %Mean solid concentration of MgO in solid 1->n steps [wt %]
    CSM.SiO2=(CSM.SiO2.*FRO+FRCSMj.SiO2)./(FR+1.0e-16);
    CSM.TiO2=(CSM.TiO2.*FRO+FRCSMj.TiO2)./(FR+1.0e-16);
    CSM.Al2O3=(CSM.Al2O3.*FRO+FRCSMj.Al2O3)./(FR+1.0e-16);
    CSM.FeO=(CSM.FeO.*FRO+FRCSMj.FeO)./(FR+1.0e-16);
    CSM.Fe2O3=(CSM.Fe2O3.*FRO+FRCSMj.Fe2O3)./(FR+1.0e-16);
    CSM.MnO=(CSM.MnO.*FRO+FRCSMj.MnO)./(FR+1.0e-16);
    CSM.MgO=(CSM.MgO.*FRO+FRCSMj.MgO)./(FR+1.0e-16);
    CSM.CaO=(CSM.CaO.*FRO+FRCSMj.CaO)./(FR+1.0e-16);
    CSM.Na2O=(CSM.Na2O.*FRO+FRCSMj.Na2O)./(FR+1.0e-16);
    CSM.K2O=(CSM.K2O.*FRO+FRCSMj.K2O)./(FR+1.0e-16);
    CSM.P2O5=(CSM.P2O5.*FRO+FRCSMj.P2O5)./(FR+1.0e-16);
    CSM.H2O=(CSM.H2O.*FRO+FRCSMj.H2O)./(FR+1.0e-16);
    %NOTE: FR=0.0 for no solid, and thus CSM=NaN, we add 1.0e-16; FRCSMj and FRCSTe are updated from MAGTFC.m.
    %MgO in solid En is 18.612 wt%

    FRO=FR;
   
    [T,FL,QM,RESE]=MAGTFC(T,FL,FKT,KL,RL,CPL,VX,VY,FRCPST,PA);
    %NOTE: MAGTFC returns T(1,:)=TW, FL(1,:)=0.0
    %MAGTFC updates T, FL, MCL.ANY, dFSLH.ANY, TCL.ANY
    %--------------------------- T-FS-CL ----------------------------------
    
    %--------------- UPDATE PARAMETERS WITH NO CONVECTION -----------------
    %Updated in MAGTFC.m: MISS, TISS, dFSLH,ANY, dFS.ANY, FS.ANY, MCS.ANY, RS.ANY, MCL.ANY, TCL.ANY, TCS.ANY, FRCSMj.ANY, KP.ANY, FRCSTe.ANY, T, FL
    
    %Update RL, CPL and MU with updated T, MCL
    [MU,CPL,RL]=VisCpRLs(T,PA,MCL,1.0-FL);
    %NOTE: MAGTFC returns T(1,:)=TW, FL(1,:)=0.0
    %NOTE: MU, CPL, RL at ghost cells have been set equal to inner hot layer, eps. top cool boundary
    
    %Update RS with new T and PA (presure range is small so that pressure effect on density is ignored!)
    RS=SolidDensitys(T,PA);%solid phase density [kg/m^3]
    %NOTE: MAGTFC returns T(1,:)=TW, FL(1,:)=0.0
    %     RST=SolidDensity(TW,PW);%TW, PW boundary solid density [kg/m^3]
    %     RS.OL(1,:)=RST.OL;%top cool boundary
    %     RS.OPX(1,:)=RST.OPX;
    %     RS.CPX(1,:)=RST.CPX;
    %     RS.PL(1,:)=RST.PL;
    %     RS.ILM(1,:)=RST.ILM;
    %Rs doesn't give OL, OPX and ILM density temporarily
    %NOTE: RS at ghost cells are set equal to inner layer since the
    %medium surrouding this cavity is not given, hence RS is not given by
    %T in ghost cells, we use Border.m in SolidDensitys.m
    
    %Update CPS.ANY
    %NOTE 1: CpS Data from Constitution of the Moon: 5. Constraints on composition,
    %density, temperature, and radius of a core, by Kuskov 1998
    %NOTE 2: T in Kelvin
    %NOTE 3: CPS at ghost cells are set equal to inner layer since the
    %medium surrouding this cavity is not given, hence CPS is not given by
    %T in ghost cells, we use Border.m.
    CPS.PL=290.90+0.0276*T-34080000.0./T.^2+5218000000.0./T.^3+29625.0./T;%[J/mol/K]
    CPS.PL=CPS.PL/0.2782073;%CaAl2Si2O8: 278.2073 g/mol [J/kg/K]
    CPS.PL=Border(CPS.PL);
    %     CPS.PL(1,:)=290.90+0.0276*TW-34080000.0/TW^2+5218000000.0/TW^3+29625.0/TW;%top boundary [J/mol/K]
    %     CPS.PL(1,:)=CPS.PL(1,:)/0.2782073;%CaAl2Si2O8: 278.2073 g/mol [J/kg/K]
    CPS.CPX=157.3+0.0000205*T-1372950.0./T.^2-1010.5./sqrt(T);%[J/mol/K]
    CPS.CPX=CPS.CPX/0.1082752;%Ca0.5Mg0.5SiO3: 108.2752 g/mol [J/kg/K]
    CPS.CPX=Border(CPS.CPX);
    %     CPS.CPX(1,:)=157.3+0.0000205*TW-1372950.0/TW^2-1010.5/sqrt(TW);%[J/mol/K]
    %     CPS.CPX(1,:)=CPS.CPX(1,:)/0.1082752;%Ca0.5Mg0.5SiO3: 108.2752 g/mol [J/kg/K]
    
    %Update KS.ANY
    KS=Conducts(T,PA);
    %NOTE: MAGTFC returns T(1,:)=TW, FL(1,:)=0.0
    %     KST=Conduct(TW,PW);%top boundary
    %     KS.OL(1,:)=KST.OL;
    %     KS.OPX(1,:)=KST.OPX;
    %     KS.CPX(1,:)=KST.CPX;
    %     KS.PL(1,:)=KST.PL;
    %     KS.ILM(1,:)=KST.ILM;
    %NOTE 1: KS is a struct
    % KS=struct('OL',zeros(NIY+2,NIX+2),...%olivine thermal conductivity [W/m/K]
    %     'OPX',zeros(NIY+2,NIX+2),...%opx thermal conductivity [W/m/K]
    %     'CPX',1.4*ones(NIY+2,NIX+2),...%cpx thermal conductivity [W/m/K]
    %     'PL',1.4*ones(NIY+2,NIX+2),...%pl thermal conductivity [W/m/K]
    %     'ILM',zeros(NIY+2,NIX+2));%ilmentite thermal conductivity [W/m/K]
    %NOTE 2: surprisingly, ilmenite has a very low thermal conductivity, see 'Thermal conductivity of titanium slags' by Heimo 2019
    %NOTE 3: KS at ghost cells are set equal to inner layer since the
    %medium surrouding this cavity is not given, hence KS is not given by
    %T in ghost cells, we use Border.m in SolidConduct.m.
    
    %Update KL
    ML.OL=ML.OL-dFSLH.OL.*RS.OL;
    ML.OPX=ML.OPX-dFSLH.OPX.*RS.OPX;
    ML.CPX=ML.CPX-dFSLH.CPX.*RS.CPX;
    ML.PL=ML.PL-dFSLH.PL.*RS.PL;
    ML.ILM=ML.ILM-dFSLH.ILM.*RS.ILM;
    %dFSLH.ANY(1,:)=0.0
    MLF.CPX=ML.CPX*100.0./(ML.OL+ML.OPX+ML.CPX+ML.PL+ML.ILM);
    MLF.PL=ML.PL*100.0./(ML.OL+ML.OPX+ML.CPX+ML.PL+ML.ILM);
    mf=(MLF.CPX/216.5504)./((MLF.CPX/216.5504)+(MLF.PL/278.2073));%Di mole fraction
    DAn=(0.36+0.4*exp(-(T-273.0)/300.0))*1.0e-6;%An thermal diffusivity [m^2/sec]
    KL=(DAn.*(1.0-mf)+0.57e-6*mf).*RL.*CPL+8.5e-11*T.^3;%molar fraction averaged thermal conductivity of liquid Klat+Krad[W/m/K]
    KL(1,:)=KL(2,:);%Psuedo top boundary
    %KL=(0.8*0.290e-6+0.2*0.358e-6)*RL.*CPL;%mass fraction averaged thermal conductivity of liquid [W/m/k]
    %NOTE 1:
    %An=0.358e-6 m^2/sec glass remelted lattice thermal diffusivity from Hofmeister 2009 (not valid)
    %Di=0.290e-6 m^2/sec glass remelted lattice thermal diffusivity from Hofmeister 2009 (not valid)
    %NOTE 2:
    %An = Table 2 of Branlund 2012 = Eq.(7) in Christopher J. Grose = (0.36+0.4*exp(-(T-273.0)/300.0))*1.0e-6
    %Di=0.57e-6 m^2/sec crystal melted to liquid lattice thermal diffusivity from Hofmeister 2008
    %NOTE: RL and CPL have been set equal at ghost cells to inner layer by
    %VisCpRLs.m, hence KL at ghost cells are equal to inner layer.
    
    %Update oxides diffusivity in liquid at new temperature and pressure
    LiqDiff(T,PA);
    %NOTE: DL at ghost cells are set equal to inner layer since the
    %medium surrouding this cavity is not given, hence DL is not given by
    %T in ghost cells
    
    %--------------- UPDATE PARAMETERS WITH NO CONVECTION -----------------
    
    %---------------------------- VX-VY-P ---------------------------------
    %New Absolute Solid Velocity [m/sec]
    %NASV(VX,VY,MU,RL);
    %This function uses new radius of solid, new FS of solid, new RS of solid, new RL of liquid, new dynamic viscosity of pure liquid and old VSXR, old VSYR
    %to give new VSXR and new VSYR, then calculates new VSX and new VSY.
    %NOTE: If solid is set static, just set VSX.ANY=0.0, VSY.ANY=0.0 and
    %FSCR1 < 0.0!
    
    %[VX,VY,PR,RESM]=VPGM(VX,VY,PR,MUO,MU,FLO,FL,RLO,RLO,RSO,RSO,dFS);
    %VPGM input parameters in order: VX=x-axis velocity, VY=y-axis velocity,
    %PR=dynamic pressure, MUO=old dynamic viscosity, MU=new dynamic viscosity,
    %FLO=old FL, FL=new FL, RLO=old RL, RL=new RL,
    %RSO=old RS, RS=new RS, dFS, RLO=old FL, RL=new RL, dFS
    %GM=Gravity Modified
    
    [VX,VY,PD,RESM]=VPBA(VX,VY,PD,MUO,MU,FLO,FL,RLB,RLB,RSB,RSB,dFS,RLO,RL);
    %[VX,VY,PR,RESM]=VPBA(VX,VY,PR,MUO,MU,FLO,FL,RLO,RL,RSO,RS,dFS,RLO,RL);
    %VPBM input parameters in order: VX=x-axis velocity, VY=y-axis velocity,
    %PR=dynamic pressure, MUO=old dynamic viscosity, MU=new dynamic viscosity,
    %FLO=old FL, FL=new FL, RLO=old RL, RL=new RL,RSO=old RS, RS=new RS,
    %dFS, RLO=old FL, RL=new RL; the last two inputs are dedicated to gravity terms
    %BA=Boussinesq Approximation
    
    MURD(:,:,n+1)=MU;%record dynamic vicosity [Pa.sec]
    MUO=MU;%Update new dynamic viscosity [Pa.sec]
    %---------------------------- VX-VY-P ---------------------------------
    
    %------------------------ RECORD VARIABLES ----------------------------
    %Record main field variables
    TRD(:,:,n+1)=T;
    FLRD(:,:,n+1)=FL;
    RLRD(:,:,n+1)=RL;
    
    FSRD.OL(:,:,n+1)=FS.OL;
    FSRD.OPX(:,:,n+1)=FS.OPX;
    FSRD.CPX(:,:,n+1)=FS.CPX;
    FSRD.PL(:,:,n+1)=FS.PL;
    FSRD.ILM(:,:,n+1)=FS.ILM;
    
    MCLRD.SiO2(:,:,n+1)=MCL.SiO2;
    MCLRD.TiO2(:,:,n+1)=MCL.TiO2;
    MCLRD.Al2O3(:,:,n+1)=MCL.Al2O3;
    MCLRD.FeO(:,:,n+1)=MCL.FeO;
    MCLRD.Fe2O3(:,:,n+1)=MCL.Fe2O3;
    MCLRD.MnO(:,:,n+1)=MCL.MnO;
    MCLRD.MgO(:,:,n+1)=MCL.MgO;
    MCLRD.CaO(:,:,n+1)=MCL.CaO;
    MCLRD.Na2O(:,:,n+1)=MCL.Na2O;
    MCLRD.K2O(:,:,n+1)=MCL.K2O;
    MCLRD.P2O5(:,:,n+1)=MCL.P2O5;
    MCLRD.H2O(:,:,n+1)=MCL.H2O;
    
    MCSRD.SiO2.OL(:,:,n+1)=MCS.SiO2.OL;
    MCSRD.SiO2.OPX(:,:,n+1)=MCS.SiO2.OPX;
    MCSRD.SiO2.CPX(:,:,n+1)=MCS.SiO2.CPX;
    MCSRD.SiO2.PL(:,:,n+1)=MCS.SiO2.PL;
    MCSRD.SiO2.ILM(:,:,n+1)=MCS.SiO2.ILM;
    
    MCSRD.TiO2.OL(:,:,n+1)=MCS.TiO2.OL;
    MCSRD.TiO2.OPX(:,:,n+1)=MCS.TiO2.OPX;
    MCSRD.TiO2.CPX(:,:,n+1)=MCS.TiO2.CPX;
    MCSRD.TiO2.PL(:,:,n+1)=MCS.TiO2.PL;
    MCSRD.TiO2.ILM(:,:,n+1)=MCS.TiO2.ILM;
    
    MCSRD.Al2O3.OL(:,:,n+1)=MCS.Al2O3.OL;
    MCSRD.Al2O3.OPX(:,:,n+1)=MCS.Al2O3.OPX;
    MCSRD.Al2O3.CPX(:,:,n+1)=MCS.Al2O3.CPX;
    MCSRD.Al2O3.PL(:,:,n+1)=MCS.Al2O3.PL;
    MCSRD.Al2O3.ILM(:,:,n+1)=MCS.Al2O3.ILM;
    
    MCSRD.FeO.OL(:,:,n+1)=MCS.FeO.OL;
    MCSRD.FeO.OPX(:,:,n+1)=MCS.FeO.OPX;
    MCSRD.FeO.CPX(:,:,n+1)=MCS.FeO.CPX;
    MCSRD.FeO.PL(:,:,n+1)=MCS.FeO.PL;
    MCSRD.FeO.ILM(:,:,n+1)=MCS.FeO.ILM;
    
    MCSRD.Fe2O3.OL(:,:,n+1)=MCS.Fe2O3.OL;
    MCSRD.Fe2O3.OPX(:,:,n+1)=MCS.Fe2O3.OPX;
    MCSRD.Fe2O3.CPX(:,:,n+1)=MCS.Fe2O3.CPX;
    MCSRD.Fe2O3.PL(:,:,n+1)=MCS.Fe2O3.PL;
    MCSRD.Fe2O3.ILM(:,:,n+1)=MCS.Fe2O3.ILM;
    
    MCSRD.MnO.OL(:,:,n+1)=MCS.MnO.OL;
    MCSRD.MnO.OPX(:,:,n+1)=MCS.MnO.OPX;
    MCSRD.MnO.CPX(:,:,n+1)=MCS.MnO.CPX;
    MCSRD.MnO.PL(:,:,n+1)=MCS.MnO.PL;
    MCSRD.MnO.ILM(:,:,n+1)=MCS.MnO.ILM;
    
    MCSRD.MgO.OL(:,:,n+1)=MCS.MgO.OL;
    MCSRD.MgO.OPX(:,:,n+1)=MCS.MgO.OPX;
    MCSRD.MgO.CPX(:,:,n+1)=MCS.MgO.CPX;
    MCSRD.MgO.PL(:,:,n+1)=MCS.MgO.PL;
    MCSRD.MgO.ILM(:,:,n+1)=MCS.MgO.ILM;
    
    MCSRD.CaO.OL(:,:,n+1)=MCS.CaO.OL;
    MCSRD.CaO.OPX(:,:,n+1)=MCS.CaO.OPX;
    MCSRD.CaO.CPX(:,:,n+1)=MCS.CaO.CPX;
    MCSRD.CaO.PL(:,:,n+1)=MCS.CaO.PL;
    MCSRD.CaO.ILM(:,:,n+1)=MCS.CaO.ILM;
    
    MCSRD.Na2O.OL(:,:,n+1)=MCS.Na2O.OL;
    MCSRD.Na2O.OPX(:,:,n+1)=MCS.Na2O.OPX;
    MCSRD.Na2O.CPX(:,:,n+1)=MCS.Na2O.CPX;
    MCSRD.Na2O.PL(:,:,n+1)=MCS.Na2O.PL;
    MCSRD.Na2O.ILM(:,:,n+1)=MCS.Na2O.ILM;
    
    MCSRD.K2O.OL(:,:,n+1)=MCS.K2O.OL;
    MCSRD.K2O.OPX(:,:,n+1)=MCS.K2O.OPX;
    MCSRD.K2O.CPX(:,:,n+1)=MCS.K2O.CPX;
    MCSRD.K2O.PL(:,:,n+1)=MCS.K2O.PL;
    MCSRD.K2O.ILM(:,:,n+1)=MCS.K2O.ILM;
    
    MCSRD.P2O5.OL(:,:,n+1)=MCS.P2O5.OL;
    MCSRD.P2O5.OPX(:,:,n+1)=MCS.P2O5.OPX;
    MCSRD.P2O5.CPX(:,:,n+1)=MCS.P2O5.CPX;
    MCSRD.P2O5.PL(:,:,n+1)=MCS.P2O5.PL;
    MCSRD.P2O5.ILM(:,:,n+1)=MCS.P2O5.ILM;
    
    MCSRD.H2O.OL(:,:,n+1)=MCS.H2O.OL;
    MCSRD.H2O.OPX(:,:,n+1)=MCS.H2O.OPX;
    MCSRD.H2O.CPX(:,:,n+1)=MCS.H2O.CPX;
    MCSRD.H2O.PL(:,:,n+1)=MCS.H2O.PL;
    MCSRD.H2O.ILM(:,:,n+1)=MCS.H2O.ILM;
    
    TCLRD.Sm(:,:,n+1)=TCL.Sm;
    TCLRD.Nd(:,:,n+1)=TCL.Nd;
    
    TCSRD.Sm.OL(:,:,n+1)=TCS.Sm.OL;
    TCSRD.Sm.OPX(:,:,n+1)=TCS.Sm.OPX;
    TCSRD.Sm.CPX(:,:,n+1)=TCS.Sm.CPX;
    TCSRD.Sm.PL(:,:,n+1)=TCS.Sm.PL;
    TCSRD.Sm.ILM(:,:,n+1)=TCS.Sm.ILM;
    
    TCSRD.Nd.OL(:,:,n+1)=TCS.Nd.OL;
    TCSRD.Nd.OPX(:,:,n+1)=TCS.Nd.OPX;
    TCSRD.Nd.CPX(:,:,n+1)=TCS.Nd.CPX;
    TCSRD.Nd.PL(:,:,n+1)=TCS.Nd.PL;
    TCSRD.Nd.ILM(:,:,n+1)=TCS.Nd.ILM;
    
    dFSRD.OL(:,:,n)=dFS.OL;
    dFSRD.OPX(:,:,n)=dFS.OPX;
    dFSRD.CPX(:,:,n)=dFS.CPX;
    dFSRD.PL(:,:,n)=dFS.PL;
    dFSRD.ILM(:,:,n)=dFS.ILM;
    
    VXRD(:,:,n+1)=VX;
    VYRD(:,:,n+1)=VY;
    
    PDRD(:,:,n+1)=PD;
    
    QMRD(n+1)=QM;
    %----------------------- RECORD VARIABLES -----------------------------
    
    %-------------------------- GRAPHICS ----------------------------------
    if(mod(n,50)==0.0)
    PlotTFC(T,FS,MCL.MgO,CSM.MgO,RL,VX,VY,PD,t(n+1),n);%VSX.PL,VSY.PL
    end
    %Current_Frame=getframe(gcf);
    %writeVideo(aviobj,Current_Frame);
    %-------------------------- GRAPHICS ----------------------------------
    
    %To find step numbers which show solidification starts (NFSS.ANY) and ends (NFSE.ANY), repectively
    %This would help if you want to know solidification location and process
    FindNSE(FSO,FS,n+1);
    
    T(1,1:NIX+2)=2.0*TW-T(2,1:NIX+2);%numerical top cool boundary
    
    fprintf('\n');
end
toc;

save('RL.mat','-mat','RL');%non-structure
save('RS.mat','-mat','RS');%structure
save('T0.mat','-mat','T');%non-structure
save('FL.mat','-mat','FL');%non-structure
save('FS.mat','-mat','FS');%structure
save('dFS.mat','-mat','dFS');%structure
save('dFSLH.mat','-mat','dFS');%structure
save('VX.mat','-mat','VX');%non-structure
save('VY.mat','-mat','VY');%non-structure
save('PD.mat','-mat','PD');%non-structure
save('FRCPS.mat','-mat','FRCPS');%structure
save('FRCSMj.mat','-mat','FRCSMj');%structure
save('FRCSTe.mat','-mat','FRCSTe');%structure
save('FK.mat','-mat','FK');%structure
save('FR.mat','-mat','FR');%non-structure
save('CPL.mat','-mat','CPL');%non-structure
save('CPS.mat','-mat','CPS');%structure
save('MCS.mat','-mat','MCS');%structure
save('TCS.mat','-mat','TCS');%structure
save('KS.mat','-mat','KS');%structure
save('KL.mat','-mat','KL');%non-structure
save('MCL.mat','-mat','MCL');%non-structure
save('TCL.mat','-mat','TCL');%non-structure
save('MU.mat','-mat','MU');%non-structure
save('MUO.mat','-mat','MUO');%non-structure
save('FLO.mat','-mat','FLO');%non-structure
save('RLO.mat','-mat','RLO');%non-structure
save('FRO.mat','-mat','FRO');%non-structure
save('MISS.mat','-mat','MISS');%structure
save('TISS.mat','-mat','TISS');%structure
save('NFSS.mat','-mat','NFSS');%structure
save('NFSE.mat','-mat','NFSE');%structure
save('CSM.mat','-mat','CSM');%structure
save('time.mat','-mat','t');%non-structure
save('DL.mat','-mat','DL');%non-structure
save('RSB.mat','-mat','RSB');%non-structure
save('RLB.mat','-mat','RLB');%non-structure
save('QM.mat','-mat','QMRD');%non-structure

%IMPORTANT NOTE: Data saved in .mat binary format will be replaced NOT APPENDED if -append option is used. While -ascii accepts APPEND option. However,
%-ascii does not support cell arrays!
time=clock;
name=[num2str(time(1)),num2str(time(2)),num2str(time(3)),num2str(time(4)),num2str(time(5))];
dirs=[pwd,'\',name];
mkdir(dirs);
cd(dirs)
save('FLRD.mat','-mat','FLRD');%structure
save('RLRD.mat','-mat','RLRD');%structure
save('MURD.mat','-mat','MURD');%structure
save('TRD.mat','-mat','TRD');%structure
save('MCLRD.mat','-mat','MCLRD');%structure
save('TCLRD.mat','-mat','TCLRD');%structure
save('dFSRD.mat','-mat','dFSRD');%structure
save('FSRD.mat','-mat','FSRD');%structure
save('VXRD.mat','-mat','VXRD');%non-structure
save('VYRD.mat','-mat','VYRD');%non-structure
save('PDRD.mat','-mat','PDRD');%non-structure
save('Time.mat','-mat','t');%non-structure

save('RL.mat','-mat','RL');%non-structure
save('RS.mat','-mat','RS');%structure
save('T0.mat','-mat','T');%non-structure
save('FL.mat','-mat','FL');%non-structure
save('FS.mat','-mat','FS');%structure
save('dFS.mat','-mat','dFS');%structure
save('dFSLH.mat','-mat','dFS');%structure
save('VX.mat','-mat','VX');%non-structure
save('VY.mat','-mat','VY');%non-structure
save('PD.mat','-mat','PD');%non-structure
save('FRCPS.mat','-mat','FRCPS');%structure
save('FRCSMj.mat','-mat','FRCSMj');%structure
save('FRCSTe.mat','-mat','FRCSTe');%structure
save('FK.mat','-mat','FK');%structure
save('FR.mat','-mat','FR');%non-structure
save('CPL.mat','-mat','CPL');%non-structure
save('CPS.mat','-mat','CPS');%structure
save('MCS.mat','-mat','MCS');%structure
save('TCS.mat','-mat','TCS');%structure
save('KS.mat','-mat','KS');%structure
save('KL.mat','-mat','KL');%non-structure
save('MCL.mat','-mat','MCL');%non-structure
save('TCL.mat','-mat','TCL');%non-structure
save('MU.mat','-mat','MU');%non-structure
save('MUO.mat','-mat','MUO');%non-structure
save('FLO.mat','-mat','FLO');%non-structure
save('RLO.mat','-mat','RLO');%non-structure
save('FRO.mat','-mat','FRO');%non-structure
save('MISS.mat','-mat','MISS');%structure
save('TISS.mat','-mat','TISS');%structure
save('NFSS.mat','-mat','NFSS');%structure
save('NFSE.mat','-mat','NFSE');%structure
save('CSM.mat','-mat','CSM');%structure
save('DL.mat','-mat','DL');%non-structure
save('RSB.mat','-mat','RSB');%non-structure
save('RLB.mat','-mat','RLB');%non-structure
save('QM.mat','-mat','QMRD');%non-structure

% fid=fopen('logs.txt');
% fprintf(fid,'Info read from  : 20213292136\n');
% fprintf(fid,'Run stopped at time : %s\n',name);
% fclose(fid);

cd ../

PlotTFC(T,FS,MCL.MgO,CSM.MgO,RL,VX,VY,PD,t(n+1),n);
OriginPlot(VX,VY,T);
%close(aviobj);


