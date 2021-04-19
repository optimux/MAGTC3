function [T,FL,QM,RESE]=MAGTFC(TTM,FLTM,FKTM,KLTM,RLTM,CPLTM,VXTM,VYTM,FRCPSTM,PATM)
%To iteratively solve energy and species conservation equations with finite
%volume method. Model, discretization and iteration can be found "A
%Continuum Model for Computer Simulation of Macrosegregations in Ingots
%during solidification" by Daming Xu 1989 and his companion papers in 1991
%Modified from TFCM.m for Xu Daming 1991
%Created on 2020-4-1
%Modified on 2020-7-1 for MAGTC3
%Modified for remelting 2021-4-9


%===== INPUT PARAMETERS ======
%TTM: temperature [K] [NIY+2,NIX+2]
%FLTM: liquid volume fraction [1] [NIY+2,NIX+2]
%FKTM: total solid thermal conductivity (total: minerals + volume fraction based) [W/m/K] [NIY+2,NIX+2]
%KLTM: liquid thermal conductivity [W/m/K] [NIY+2,NIX+2]
%RLTM: liquid density [kg/m^3] [NIY+2,NIX+2]
%CPLTM: liquid specific heat capacity [J/kg/K] [NIY+2,NIX+2]
%VSXTM: absolute solid velocity [m/sec] [NIY+2,NIX+1] structure
%VSYTM: absolute solid velocity [m/sec] [NIY+1,NIX+2] structure
%FRCPSTM: temporary FRCPS=FS*RS*CPS [J/m^3/K] [NIY+2,NIX+2] structure

global NIX
global NIY
global dx
global dy
global CE
global dtb
global TE
global FSE
global FSELOG
global MCL
global MCS
global TCL
global TCS
global TW
global FRCSMj
global FRCSTe
% global Majors
% global Minors
global RS
%RS has never been updated to RSN in this function
global RSN
%RSN has been updated by itself (i.e., RSN(iters_n+1)=RSN(iters_n)) in this function iteratively
global FS
global dFS
global CPS
%CPS has never been updated in this function
% global RFVX
% global RFVY
global HS
global dFSLH
global MCSN
global RLE
global CPLE
global RSE
global CPSE
global MISS
global TISS
global KP
global DL
global KST
global CL0
global TL0

err0dFSLH=1.0e-9;%controlled accuracy for dFSLH
err0TTM=1.0e-9;%controlled accuracy for T
err0CLTM=1.0e-9;%controlled accuracy for CL
NdFS=500;%max iterations allowed of dFS, T, CL
phi=0.5;%integration coefficient [0,1]
acc=0.5;%accelerator factor for dFS, TN, CLN,FLN

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================

%################################################################### ENERGY BALANCE #########################################################################
%% ================= DIFFUSION HEAT FLUX =========================

%------------- thermal conductivity on finite volume faces ----------------
%....................... all faces, x-axis [W/m/K] ........................
KX=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    for j=1:NIY+2
        KX(j,i)=0.5*(FKTM(j,i)+FKTM(j,i+1)+FLTM(j,i)*KLTM(j,i)+FLTM(j,i+1)*KLTM(j,i+1));%[(FS*KS)+(FL*KL)]|(j+/-0.5,k)
    end
end
%........................all faces, x-axis [W/m/K] ........................

%........................all faces, y-axis [W/m/K] ........................
KY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        KY(j,i)=0.5*(FKTM(j,i)+FKTM(j+1,i)+FLTM(j,i)*KLTM(j,i)+FLTM(j+1,i)*KLTM(j+1,i));%[(FS*KS)+(FL*KL)]|(j,k+/-0.5)
    end
end
%KY(1,:)= top real boundary at T_radiation, constant composition and varied temperature, so its thermal conductivity changes with time
%NOTE: first row is complete solid, it serves as environment and KY(1,:) should be K of complete solid K(1,:). However, we use (K(1,:)+K(2,:))/2 as smoothed
%average.
%........................all faces, y-axis [W/m/K] ........................

%------------- thermal conductivity on finite volume faces ----------------

%----------------------------- T Gradient ---------------------------------

DTX=zeros(NIY+2,NIX+1);%T gradient between finite volumes in x-axis
for i=2:NIX
    for j=1:NIY+2
        DTX(j,i)=(TTM(j,i+1)-TTM(j,i))/(0.5*dx(i-1)+0.5*dx(i));%[K/m]
    end
end

%Recheck boundary condition explicitly
for i=1:NIY+2
    %TTM(1:NIY+2,1)==TTM(1:NIY+2,2) --> DTX(1:NIY+2,1)=0.0
    DTX(i,1)=(TTM(i,2)-TTM(i,1))/(0.5*dx(1)+0.5*dx(1));%left insulated boundary
    
    %TTM(1:NIY+2,NIX+2)==TTM(1:NIY+2,NIX+1) --> DTX(1:NIY+2,NIX+1)=0.0
    DTX(i,NIX+1)=(TTM(i,NIX+2)-TTM(i,NIX+1))/(0.5*dx(NIX)+0.5*dx(NIX));%right insulated boundary
end

DTY=zeros(NIY+1,NIX+2);%T gradient between finite volumes in y-axis
for i=1:NIX+2
    for j=2:NIY
        DTY(j,i)=(TTM(j+1,i)-TTM(j,i))/(0.5*dy(j-1)+0.5*dy(j));%[K/m]
    end
end

%Recheck boundary condition explicitly
for i=1:NIX+2
    %TTM(1,1:NIX+2)!=TTM(2,1:NIX+2) --> DTY(1,1:NIX+2)!=0.0
    DTY(1,i)=(TTM(2,i)-TTM(1,i))/(0.5*dy(1)+0.5*dy(1));%top cool boundary
    
    %TTM(NIY+2,1:NIX+2)==TTM(NIY+1,1:NIX+2) --> DTY(NIY+1,1:NIX+1)=0
    DTY(NIY+1,i)=(TTM(NIY+2,i)-TTM(NIY+1,i))/(0.5*dy(NIY)+0.5*dy(NIY));%bottom insulated boundary
    %In Spera 1997, bottom temperature is set as constant.
end

%----------------------------- T gradient ---------------------------------

%----------------------------- Heat flux ----------------------------------

%heat flux on all x faces
TQDX=zeros(NIY+2,NIX+1);%T=T, Q=flux, D=Diffusion, X=x-axis
for i=1:NIX+1
    for j=1:NIY+2
        TQDX(j,i)=KX(j,i)*DTX(j,i);%[W/m^2]
    end
end

%heat flux on all y faces
TQDY=zeros(NIY+1,NIX+2);%T=T, Q=flux, D=Diffusion, Y=y-axis
for i=1:NIX+2
    for j=1:NIY+1
        TQDY(j,i)=KY(j,i)*DTY(j,i);%[W/m^2]
    end
end

%total heat flux of each finite volume
TQDT=zeros(NIY,NIX);%T=T, Q=flux, D=Diffusion, T=total
for i=1:NIX
    for j=1:NIY
        TQDT(j,i)=-(TQDX(j+1,i)-TQDX(j+1,i+1))*dtb/dx(i)-(TQDY(j,i+1)-TQDY(j+1,i+1))*dtb/dy(j);%[J/m^3]
    end
end

%----------------------------- Heat flux ----------------------------------

%% =================== INTERNAL HEAT STORAGE =====================
IHS=zeros(NIY+2,NIX+2);%old T related Internal Heat Storage

for i=1:NIX+2
    for j=1:NIY+2
        IHS(j,i)=FRCPSTM(j,i)*TTM(j,i)+FLTM(j,i)*RLTM(j,i)*CPLTM(j,i)*TTM(j,i);%[J/m^3]
    end
end

%% ================ LIQUID CONVECTION HEAT FLUX ==================

%------------------------ Part One: RFVX RFVY -----------------------------
%                          mass flux density

%x-axis velocity
%        +-------+
%        |       |
%    --> |   *   |  -->
%        |       |
%        +-------+
%NOTE: in fact, vx at top and bottom faces should be estimated,
%but they can be interpolated; so only vx normal to vertical faces are
%calculated!

RFVX=zeros(NIY+2,NIX+1);%RL*FL*VX, all VX at volume faces
for i=2:NIX
    %VX(1:NIY+2,1)=0.0 --> RFVX(1:NIY+2,1)=0.0 (left impermeable boundary)
    %VX(1:NIY+2,NIX+1)=0.0 --> RFVX(1:NIY+2,NIX+1)=0.0 (right impermeable boundary)
    for j=1:NIY+2
        %VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> RFVX(1,1:NIX+1)=-RFVX(2,1:NIX+1) (top NO SLIP boundary)
        %VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1) --> RFVX(NIY+2,1:NIX+1)=-RFVX(NIY+1,1;NIX+1) (bottom NO SLIP boundary)
        
        %VX(1,1:NIX+1)=VX(2,1:NIX+1) --> RFVX(1,1:NIX+1)=RFVX(2,1:NIX+1) (top FREE boundary)
        %VX(NIY+2,1:NIX+1)=VX(NIY+1,1:NIX+1) --> RFVX(NIY+2,1:NIX+1)=RFVX(NIY+1,1;NIX+1) (bottom FREE boundary)
        RFVX(j,i)=VXTM(j,i)*0.5*(RLTM(j,i)+RLTM(j,i+1))*0.5*(FLTM(j,i)+FLTM(j,i+1));%[kg/m^2/sec]
    end
end

%Recheck boundary condition explicitly
for i=1:NIX+1
    RFVX(1,i)=-RFVX(2,i);%top NO SLIP boundary
    RFVX(NIY+2,i)=RFVX(NIY+1,i);%bottom FREE boundary
end

% for i=1:NIX+1
% RFVX(1,i)=RFVX(2,i);%top FREE boundary
% RFVX(NIY+2,i)=RFVX(NIY+1,i);%bottom FREE boundary
% end

%y-axis velocity
%            ^
%            |
%        +-------+
%        |       |
%        |   *   |
%        |       |
%        +-------+
%            ^
%            |
%NOTE: in fact, vy at left and right faces should be estimated,
%but they can be interpolated; so only vy normal to horizontal faces are
%calculated!

RFVY=zeros(NIY+1,NIX+2);%RL*FL*VY, all VY faces
for i=1:NIX+2
    %VY(1:NIY+1,1)=VY(1:NIY+1,2) --> RFVY(1:NIY+1,1)=RFVY(1:NIY+1,2) (left FREE VY)
    %VY(1:NIY+1,NIX+2)=VY(1:NIY+1,NIX+1) --> RFVY(1:NIY+1,NIX+2)=RFVY(1:NIY+1,NIX+1) (right FREE VY)
    
    %VY(1:NIY+1,1)=-VY(1:NIY+1,2) --> RFVY(1:NIY+1,1)=-RFVY(1:NIY+1,2) (left NO SLIP VY)
    %VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1) --> RFVY(1:NIY+1,NIX+2)=-RFVY(1:NIY+1,NIX+1) (right NO SLIP VY)
    
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> RFVY(1,1:NIX+2)=0.0 (top impermeable boundary)
        %VY(NIY+1,1:NIX+2)=0.0 --> RFVY(NIY+1,1:NIX+2)=0.0 (bottom impermeable boundary)
        RFVY(j,i)=VYTM(j,i)*0.5*(RLTM(j,i)+RLTM(j+1,i))*0.5*(FLTM(j,i)+FLTM(j+1,i));%[kg/m^2/sec]
    end
    
end

%Recheck boundary condition explicitly
for i=1:NIY+1
    RFVY(i,1)=RFVY(i,2);%left FREE boundary
    RFVY(i,NIX+2)=RFVY(i,NIX+1);%right FREE boundary
end

% for i=1:NIY+1
%     RFVY(i,1)=-RFVY(i,2);%left NO SLIP boundary
%     RFVY(i,NIX+2)=-RFVY(i,NIX+1);%right NO SLIP boundary
% end
%------------------------ Part One: RFVX RFVY -----------------------------

%----------------------- Part Two: TRFVX TRFVY ----------------------------
TRFVX=zeros(NIY+2,NIX+1);
for i=2:NIX
    %VX(1:NIY+2,1)=0.0 --> TRFVX(1:NIY+2,1)=0.0 (left impermeable boundary)
    %VX(1:NIY+2,NIX+1)=0.0 --> TRFVX(1:NIY+2,NIX+1)=0.0 (right impermeable boundary)
    for j=1:NIY+2
        TRFVX(j,i)=TTM(j,i)*max(RFVX(j,i),0.0)+TTM(j,i+1)*min(RFVX(j,i),0.0);%[kg.K/m^2/sec]
    end
    %TRFVX(1,:) and TRFVX(NIY+2,:) is useless, so its value doesn't matter
end

TRFVY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> TRFVY(1,1:NIX+2)=0.0 (top impermeable boundary )
        %VY(NIY+1,1:NIX+2)=0.0 --> TRFVY(NIY+1,1:NIX+2)=0.0 (bottom impermeable boundary)
        TRFVY(j,i)=TTM(j,i)*max(RFVY(j,i),0.0)+TTM(j+1,i)*min(RFVY(j,i),0.0);%[kg.K/m^2/sec]
    end
    %TRFVY(1,:) and TRFVY(NIY+1,:) is useless, so its value doesn't matter
end

%Recheck boundary condition explicitly
for i=1:NIY+1
    TRFVY(i,1)=TRFVY(i,2);%left FREE boundary --> Some heat transfered up and down by mass flow at left boundary
    TRFVY(i,NIX+2)=TRFVY(i,NIX+1);%right FREE boundary --> Some heat transfered up and down by mass flow at right boundary
end

% for i=1:NIY+1
%     TRFVY(i,1)=-TRFVY(i,2);%left NO SLIP boundary --> No heat transfered up and down by mass flow at left boundary
%     TRFVY(i,NIX+2)=-TRFVY(i,NIX+1);%right NO SLIP boundary --> No heat transfered up and down by mass flow at right boundary
% end
%----------------------- Part Two: TRFVX TRFVY ----------------------------

%----------------------- Part Three: Heat flux ----------------------------
TQVX=zeros(NIY+2,NIX);%T=T, Q=flux, VX=VX
for i=1:NIX
    for j=1:NIY+2
        TQVX(j,i)=-dtb*CPLTM(j,i+1)*(TRFVX(j,i+1)-TRFVX(j,i))/dx(i);%[J/m^3]
    end
    %TQVX(1,:) and TQVX(NIY+2,:) is useless, its value doesn't matter
end

TQVY=zeros(NIY,NIX+2);%T=T, Q=flux, VY=VY
for i=1:NIX+2
    for j=1:NIY
        TQVY(j,i)=-dtb*CPLTM(j+1,i)*(TRFVY(j+1,i)-TRFVY(j,i))/dy(j);%[J/m^3]
    end
    %TQVY(:,1) and TQVY(:,NIX+2) is useless, its value doesn't matter
end

TQVT=zeros(NIY,NIX);%T=T, Q=flux, V=velocity, T=total
TQVT(1:NIY,1:NIX)=TQVX(2:NIY+1,1:NIX)+TQVY(1:NIY,2:NIX+1);%[J/m^3]

%------------------------ Part Three: Heat flux ---------------------------

%% ================ SOILD CONVECTION HEAT FLUX ===================

%------------------------- Part Three:  Heat flux --------------------------

%################################################################## ENERGY BALANCE #######################################################################


%################################################################## SPECIES BALANCE ######################################################################

CLQDT=struct('SiO2',zeros(NIY,NIX),...
    'TiO2',zeros(NIY,NIX),...
    'Al2O3',zeros(NIY,NIX),...
    'FeO',zeros(NIY,NIX),...
    'Fe2O3',zeros(NIY,NIX),...
    'MnO',zeros(NIY,NIX),...
    'MgO',zeros(NIY,NIX),...
    'CaO',zeros(NIY,NIX),...
    'Na2O',zeros(NIY,NIX),...
    'K2O',zeros(NIY,NIX),...
    'P2O5',zeros(NIY,NIX),...
    'H2O',zeros(NIY,NIX),...
    'Sm',zeros(NIY,NIX),...
    'Nd',zeros(NIY,NIX));%CL=CL, Q=flux, D=diffusion, T=total [kg/m^3]

% ISS=struct('SiO2',zeros(NIY+2,NIX+2),...
%     'TiO2',zeros(NIY+2,NIX+2),...
%     'Al2O3',zeros(NIY+2,NIX+2),...
%     'FeO',zeros(NIY+2,NIX+2),...
%     'Fe2O3',zeros(NIY+2,NIX+2),...
%     'MnO',zeros(NIY+2,NIX+2),...
%     'MgO',zeros(NIY+2,NIX+2),...
%     'CaO',zeros(NIY+2,NIX+2),...
%     'Na2O',zeros(NIY+2,NIX+2),...
%     'K2O',zeros(NIY+2,NIX+2),...
%     'P2O5',zeros(NIY+2,NIX+2),...
%     'H2O',zeros(NIY+2,NIX+2),...
%     'Sm',zeros(NIY+2,NIX+2),...
%     'Nd',zeros(NIY+2,NIX+2));%Internal Species Storage [kg/m^3]

CLQVT=struct('SiO2',zeros(NIY,NIX),...
    'TiO2',zeros(NIY,NIX),...
    'Al2O3',zeros(NIY,NIX),...
    'FeO',zeros(NIY,NIX),...
    'Fe2O3',zeros(NIY,NIX),...
    'MnO',zeros(NIY,NIX),...
    'MgO',zeros(NIY,NIX),...
    'CaO',zeros(NIY,NIX),...
    'Na2O',zeros(NIY,NIX),...
    'K2O',zeros(NIY,NIX),...
    'P2O5',zeros(NIY,NIX),...
    'H2O',zeros(NIY,NIX),...
    'Sm',zeros(NIY,NIX),...
    'Nd',zeros(NIY,NIX));%CL=CL, Q=flux, V=liquid velocity, T=total [kg/m^3]

% CSQVT=struct('SiO2',zeros(NIY,NIX),...
%     'TiO2',zeros(NIY,NIX),...
%     'Al2O3',zeros(NIY,NIX),...
%     'FeO',zeros(NIY,NIX),...
%     'Fe2O3',zeros(NIY,NIX),...
%     'MnO',zeros(NIY,NIX),...
%     'MgO',zeros(NIY,NIX),...
%     'CaO',zeros(NIY,NIX),...
%     'Na2O',zeros(NIY,NIX),...
%     'K2O',zeros(NIY,NIX),...
%     'P2O5',zeros(NIY,NIX),...
%     'H2O',zeros(NIY,NIX),...
%     'Sm',zeros(NIY,NIX),...
%     'Nd',zeros(NIY,NIX));%CS=CS, Q=flux, V=solid velocity, T=total

[CLQDT.SiO2,CLQVT.SiO2]=Species(MCL.SiO2,DL.SiO2,FLTM,RLTM,RFVX,RFVY);
[CLQDT.TiO2,CLQVT.TiO2]=Species(MCL.TiO2,DL.TiO2,FLTM,RLTM,RFVX,RFVY);
[CLQDT.Al2O3,CLQVT.Al2O3]=Species(MCL.Al2O3,DL.Al2O3,FLTM,RLTM,RFVX,RFVY);
[CLQDT.FeO,CLQVT.FeO]=Species(MCL.FeO,DL.FeO,FLTM,RLTM,RFVX,RFVY);
[CLQDT.Fe2O3,CLQVT.Fe2O3]=Species(MCL.Fe2O3,DL.Fe2O3,FLTM,RLTM,RFVX,RFVY);
[CLQDT.MnO,CLQVT.MnO]=Species(MCL.MnO,DL.MnO,FLTM,RLTM,RFVX,RFVY);
[CLQDT.MgO,CLQVT.MgO]=Species(MCL.MgO,DL.MgO,FLTM,RLTM,RFVX,RFVY);
[CLQDT.CaO,CLQVT.CaO]=Species(MCL.CaO,DL.CaO,FLTM,RLTM,RFVX,RFVY);
[CLQDT.Na2O,CLQVT.Na2O]=Species(MCL.Na2O,DL.Na2O,FLTM,RLTM,RFVX,RFVY);
[CLQDT.K2O,CLQVT.K2O]=Species(MCL.K2O,DL.K2O,FLTM,RLTM,RFVX,RFVY);
[CLQDT.P2O5,CLQVT.P2O5]=Species(MCL.P2O5,DL.P2O5,FLTM,RLTM,RFVX,RFVY);
[CLQDT.H2O,CLQVT.H2O]=Species(MCL.H2O,DL.H2O,FLTM,RLTM,RFVX,RFVY);
[CLQDT.Sm,CLQVT.Sm]=Species(TCL.Sm,DL.Sm,FLTM,RLTM,RFVX,RFVY);
[CLQDT.Nd,CLQVT.Nd]=Species(TCL.Nd,DL.Nd,FLTM,RLTM,RFVX,RFVY);

%################################################################## SPECIES BALANCE ######################################################################

%################################################################ T-FS-CL ITERATION ########################################################################

TFN=zeros(NIY+2,NIX+2);%T Factor of Next step
TN=zeros(NIY+2,NIX+2);%T of Next step
TLTM=zeros(NIY+2,NIX+2);%Liquidus of old step
iters=zeros(NIY+2,NIX+2);
TTMBE=TTM;%old T before eutectic for eutectic T calculation

RSCPSdFSTM=zeros(NIY+2,NIX+2);
dFSLHN=struct('OL',zeros(NIY+2,NIX+2),...%New time step olivine volume fraction increment related to latent heat [1]
    'OPX',zeros(NIY+2,NIX+2),...%New time step opx volume fraction increment related to latent heat [1]
    'CPX',zeros(NIY+2,NIX+2),...%New time step cpx volume fraction increment related to latent heat [1]
    'PL',zeros(NIY+2,NIX+2),...%New time step pl volume fraction increment related to latent heat [1]
    'ILM',zeros(NIY+2,NIX+2));%New time step ilmentite volume fraction increment related to latent heat [1]
%dFSLHN.ANY(1,:) must be 0.0 to match real top cold boundary since RSCPSdFSTM, dFSLHT, TFN, LH.ANY and LHT depend on dFSLHN.ANY
dFSLHT=zeros(NIY+2,NIX+2);%total solid volume fraction increment related to latent heat
LHT=zeros(NIY+2,NIX+2);%total latent heat release in action [J/m^3]

LH=struct('OL',zeros(NIY+2,NIX+2),...%olivine latent heat released [J/m^3]
    'OPX',zeros(NIY+2,NIX+2),...%opx latent heat released [J/m^3]
    'CPX',zeros(NIY+2,NIX+2),...%cpx latent heat released [J/m^3]
    'PL',zeros(NIY+2,NIX+2),...%pl latent heat released [J/m^3]
    'ILM',zeros(NIY+2,NIX+2));%ilmentite latent heat released [J/m^3]

%TFN: T Factor of Next step
for i=1:NIX+2
    for j=1:NIY+2
        RSCPSdFSTM(j,i)=RS.OL(j,i)*CPS.OL(j,i)*dFSLHN.OL(j,i)+...
            RS.OPX(j,i)*CPS.OPX(j,i)*dFSLHN.OPX(j,i)+...
            RS.CPX(j,i)*CPS.CPX(j,i)*dFSLHN.CPX(j,i)+...
            RS.PL(j,i)*CPS.PL(j,i)*dFSLHN.PL(j,i)+...
            RS.ILM(j,i)*CPS.ILM(j,i)*dFSLHN.ILM(j,i);%[J/m^3/K]
        
        dFSLHT(j,i)=dFSLHN.OL(j,i)+dFSLHN.OPX(j,i)+dFSLHN.CPX(j,i)+dFSLHN.PL(j,i)+dFSLHN.ILM(j,i);%[1]
        TFN(j,i)=FRCPSTM(j,i)+FLTM(j,i)*RLTM(j,i)*CPLTM(j,i)+RSCPSdFSTM(j,i)-RLTM(j,i)*CPLTM(j,i)*dFSLHT(j,i);%[J/m^3/K]
        
        LH.OL(j,i)=dFSLHN.OL(j,i)*RS.OL(j,i)*HS.OL(j,i);%[J/m^3]
        LH.OPX(j,i)=dFSLHN.OPX(j,i)*RS.OPX(j,i)*HS.OPX(j,i);
        LH.CPX(j,i)=dFSLHN.CPX(j,i)*RS.CPX(j,i)*HS.CPX(j,i);
        LH.PL(j,i)=dFSLHN.PL(j,i)*RS.PL(j,i)*HS.PL(j,i);
        LH.ILM(j,i)=dFSLHN.ILM(j,i)*RS.ILM(j,i)*HS.ILM(j,i);
        
        LHT(j,i)=LH.OL(j,i)+LH.OPX(j,i)+LH.CPX(j,i)+LH.PL(j,i)+LH.ILM(j,i);%[J/m^3]
        
    end
end
%RSCPSdFSTM: RS*CPS at n step, as approxiamation of RS*CPS at n+1 step
%RLTM, CPLTM: RS*CPL at n step, as approxiamation of RS*CPL at n+1 step

%----------------------------- STEP ONE -----------------------------------
%energy loss due to diffusion and convection [J/m^3]
dQ=zeros(NIY+2,NIX+2);
%Let the initial value of dFSLH.ANY equal zero and calculate the approximation of TN.
for i=2:NIX+1
    for j=2:NIY+1
        dQ(j,i)=TQDT(j-1,i-1)+TQVT(j-1,i-1);%[J/m^3]
        TN(j,i)=(dQ(j,i)+IHS(j,i)+LHT(j,i))/TFN(j,i);%[J/m^3 divided by J/m^3/K == K]
    end
end

TN(1,2:NIX+1)=2.0*TW-TN(2,2:NIX+1);%numerical top cool boundary
TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%bottom insulated boundary
TN(1:NIY+2,1)=TN(1:NIY+2,2);%left insulated boundary
TN(1:NIY+2,NIX+2)=TN(1:NIY+2,NIX+1);%right insulted boundary
TTM=TN;
TTMAE=TN;%New T after eutectic for eutectic T calculation
%----------------------------- STEP ONE -----------------------------------

%Temporary CL of major elements in liquid
CLTM=struct('SiO2',zeros(NIY+2,NIX+2),...%SiO2 in liquid [1]
    'TiO2',zeros(NIY+2,NIX+2),...%TiO2 in liquid [1]
    'Al2O3',zeros(NIY+2,NIX+2),...%Al2O3 in liquid [1]
    'FeO',zeros(NIY+2,NIX+2),...%FeO in liquid [1]
    'Fe2O3',zeros(NIY+2,NIX+2),...%Fe2O3 in liquid [1]
    'MnO',zeros(NIY+2,NIX+2),...%MnO in liquid [1]
    'MgO',zeros(NIY+2,NIX+2),...%MgO in liquid [1]
    'CaO',zeros(NIY+2,NIX+2),...%CaO in liquid [1]
    'Na2O',zeros(NIY+2,NIX+2),...%Na2O in liquid [1]
    'K2O',zeros(NIY+2,NIX+2),...%K2O in liquid [1]
    'P2O5',zeros(NIY+2,NIX+2),...%P2O5 in liquid [1]
    'H2O',zeros(NIY+2,NIX+2),...%H2O in liquid [1]
    'Sm',zeros(NIY+2,NIX+2),...%Sm in liquid [1]
    'Nd',zeros(NIY+2,NIX+2));%Nd in liquid [1]

CLTM=MCL;

%Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1; not wt%]
CLCM=MCL.MgO;
%NOTE: CLCM will be used to determine solidification stages

TLTM=-1.3351*(CLCM*100.0).^2+55.878*CLCM*100.0+813.9991+273.15;%old liquidus to determine superliquidus or subliquidus [K]

%logical variable for while loop
errdFS=ones(NIY+2,NIX+2);%initial error of dFS, to triger while loop
errTTM=ones(NIY+2,NIX+2);
errCLTM=ones(NIY+2,NIX+2);
err=zeros(NIY+2,NIX+2);
FSTTM=zeros(NIY+2,NIX+2);%old step total solid volume fraction of each cell

for i=1:NIX+2
    for j=1:NIY+2
        err(j,i)=logical((err0dFSLH<=errdFS(j,i))||(err0TTM<=errTTM(j,i))||(err0CLTM<=errCLTM(j,i)));
        FSTTM(j,i)=FS.OL(j,i)+FS.OPX(j,i)+FS.CPX(j,i)+FS.PL(j,i)+FS.ILM(j,i);
    end
end

%CL Factor of Next step
CLFN=struct('SiO2',zeros(NIY,NIX),...
    'TiO2',zeros(NIY,NIX),...
    'Al2O3',zeros(NIY,NIX),...
    'FeO',zeros(NIY,NIX),...
    'Fe2O3',zeros(NIY,NIX),...
    'MnO',zeros(NIY,NIX),...
    'MgO',zeros(NIY,NIX),...
    'CaO',zeros(NIY,NIX),...
    'Na2O',zeros(NIY,NIX),...
    'K2O',zeros(NIY,NIX),...
    'P2O5',zeros(NIY,NIX),...
    'H2O',zeros(NIY,NIX),...
    'Sm',zeros(NIY,NIX),...
    'Nd',zeros(NIY,NIX));

%CS factor of Next step
CLFNS=struct('SiO2',zeros(NIY,NIX),...
    'TiO2',zeros(NIY,NIX),...
    'Al2O3',zeros(NIY,NIX),...
    'FeO',zeros(NIY,NIX),...
    'Fe2O3',zeros(NIY,NIX),...
    'MnO',zeros(NIY,NIX),...
    'MgO',zeros(NIY,NIX),...
    'CaO',zeros(NIY,NIX),...
    'Na2O',zeros(NIY,NIX),...
    'K2O',zeros(NIY,NIX),...
    'P2O5',zeros(NIY,NIX),...
    'H2O',zeros(NIY,NIX),...
    'Sm',zeros(NIY,NIX),...
    'Nd',zeros(NIY,NIX));


%Major CL of Next step
CLN=struct('SiO2',zeros(NIY+2,NIX+2),...%SiO2 in liquid [1]
    'TiO2',zeros(NIY+2,NIX+2),...%TiO2 in liquid [1]
    'Al2O3',zeros(NIY+2,NIX+2),...%Al2O3 in liquid [1]
    'FeO',zeros(NIY+2,NIX+2),...%FeO in liquid [1]
    'Fe2O3',zeros(NIY+2,NIX+2),...%Fe2O3 in liquid [1]
    'MnO',zeros(NIY+2,NIX+2),...%MnO in liquid [1]
    'MgO',zeros(NIY+2,NIX+2),...%MgO in liquid [1]
    'CaO',zeros(NIY+2,NIX+2),...%CaO in liquid [1]
    'Na2O',zeros(NIY+2,NIX+2),...%Na2O in liquid [1]
    'K2O',zeros(NIY+2,NIX+2),...%K2O in liquid [1]
    'P2O5',zeros(NIY+2,NIX+2),...%P2O5 in liquid [1]
    'H2O',zeros(NIY+2,NIX+2),...%H2O in liquid [1]
    'Sm',zeros(NIY+2,NIX+2),...%Sm in liquid [1]
    'Nd',zeros(NIY+2,NIX+2));%Nd in liquid [1]

%Temporary MCL for VisCpRL.m
MCL0=struct('SiO2',0.0,...%SiO2 in liquid [1]
    'TiO2',0.0,...%TiO2 in liquid [1]
    'Al2O3',0.0,...%Al2O3 in liquid [1]
    'FeO',0.0,...%FeO in liquid [1]
    'Fe2O3',0.0,...%Fe2O3 in liquid [1]
    'MnO',0.0,...%MnO in liquid [1]
    'MgO',0.0,...%MgO in liquid [1]
    'CaO',0.0,...%CaO in liquid [1]
    'Na2O',0.0,...%Na2O in liquid [1]
    'K2O',0.0,...%K2O in liquid [1]
    'P2O5',0.0,...%P2O5 in liquid [1]
    'H2O',0.0);%H2O in liquid [1]

FLN=zeros(NIY+2,NIX+2);%FL of Next step
RLN=zeros(NIY+2,NIX+2);%RL of Next step
TLN=zeros(NIY+2,NIX+2);%Liquidus of Next step
TEC=zeros(NIY+2,NIX+2);%New T for Energy Check

dFSLHTM=struct('OL',zeros(NIY+2,NIX+2),...%Temporary olivine volume fraction increment related to latent heat [1]
    'OPX',zeros(NIY+2,NIX+2),...%Temporary opx volume fraction increment related to latent heat [1]
    'CPX',zeros(NIY+2,NIX+2),...%Temporary cpx volume fraction increment related to latent heat [1]
    'PL',zeros(NIY+2,NIX+2),...%Temporary pl volume fraction increment related to latent heat [1]
    'ILM',zeros(NIY+2,NIX+2));%Temporary ilmentite volume fraction increment related to latent heat [1]

STLL=0.0;%Sensible heat Transform into Latent heat, temporarily used as Left hand side
STLR=0.0;%Sensible heat Transform into Latent heat, temporarily used as Right hand side

caseID=3*ones(NIY,NIX);%case ID
%dFSSEN=zeros(NIY+2,NIX+2);%used in energy check for sensible heat

%RSE=1.0/(0.599/RSE.CPX+0.401/RSE.PL);%eutectic mean solid density [kg/m^3]
%RSEHSE=0.599*RSE*HS.CPX+0.401*RSE*HS.PL;%total RSE*HSE

%Prapare for latent heat for energy check [J/m^3]
LAT=zeros(NIY,NIX);
%Prapare for sensible heat for energy check [J/m^3]
SEN=zeros(NIY,NIX);%sensible heat

modifymarker=0;

%New RS of Single point
RSNS=struct('OL',0.0,...%olivine density [kg/m^3]
    'OPX',0.0,...%opx density [kg/m^3]
    'CPX',0.0,...%cpx density [kg/m^3]
    'PL',0.0,...%pl density [kg/m^3]
    'ILM',0.0);%ilmentite density [kg/m^3]

dT=zeros(NIY+2,NIX+2);%CPX remelting to dT [K]
%------------------------ T-FS-CL MAIN LOOP ---------------------
for i=1:NIX
    % i;
    for j=1:NIY
        % j;
        QS=FS.CPX(j+1,i+1)*RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1);%energy necessary to remelt CPX [J/m^3]
        %Prepare Temporary MCL for VisCpRL.m
        MCL0.SiO2=MCL.SiO2(j+1,i+1);
        MCL0.TiO2=MCL.TiO2(j+1,i+1);
        MCL0.Al2O3=MCL.Al2O3(j+1,i+1);
        MCL0.FeO=MCL.FeO(j+1,i+1);
        MCL0.Fe2O3=MCL.Fe2O3(j+1,i+1);
        MCL0.MnO=MCL.MnO(j+1,i+1);
        MCL0.MgO=MCL.MgO(j+1,i+1);
        MCL0.CaO=MCL.CaO(j+1,i+1);
        MCL0.Na2O=MCL.Na2O(j+1,i+1);
        MCL0.K2O=MCL.K2O(j+1,i+1);
        MCL0.P2O5=MCL.P2O5(j+1,i+1);
        MCL0.H2O=MCL.H2O(j+1,i+1);
        
        if(TN(j+1,i+1)<TE)
            dFSLHPL=(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1))*(TE-TN(j+1,i+1))/...%total heat
                (0.599*RSE.PL*HS.CPX(j+1,i+1)/0.401+RSE.PL*HS.PL(j+1,i+1)+...%latent heat
                -0.599*RSE.PL*CPSE.CPX*TE/0.401-RSE.PL*CPSE.PL*TE+(1.0+0.599*RSE.PL/(0.401*RSE.CPX))*RLTM(j+1,i+1)*CPLTM(j+1,i+1)*TE);
            %NOTE: In solid property change part, we use TE because TE>TN. TE>TN means dFSLHPL becomes smaller, and if smaller dFSLHPL+FSTTM>1.0, then
            %dFSLHPL(TN)+FSTTM must be larger, that is, we use a smaller dFSLHPL+FSTTM to decide caseID
            
            %dFSLHPL=(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1))*(TE-TN(j+1,i+1))/(0.599*RSE.PL*HS.CPX(j+1,i+1)/0.401+RSE.PL*HS.PL(j+1,i+1));%RLE==RLTM, CPLE==CPLTM
            dFSLHCPX=0.599*dFSLHPL*RSE.PL/(0.401*RSE.CPX);
            dFSLHCAL=dFSLHPL+dFSLHCPX;
            
            %ORIGINALLY AS: dFSCAL=(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1))*(TE-TN(j+1,i+1))/(RSEHSE-(RSCPSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1))*TE);%Eq.(14)
            %This equation is a first guess of dFSLH in SITUATION 4&5 where T(n)=TE and T(n+1)_first_guess=TN <TE
            %(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1))*(TE-TN(j+1,i+1)) is minus total energy loss of system;
            %RSE*HSE is heat source while RSCPSTM-RLCPLTM is heat sink since the heat storage capacity of solid
            %is better than that of liquid, meaning solid keeps more energy
        else
            dFSLHPL=0.0;
            dFSLHCPX=0.0;
            dFSLHCAL=0.0;
        end
        
        %IMPORTANT NOTE: Compare T, FS and CLCM with TE, FSE and CE to determine Flags of condition
        
        %SITUATION 3: MUSHY --> EUTECTIC
        %0<FSTTM(j+1,i+1)<FSE, TTM9j+1,i+1)==TLTM(j+1,i+1)>TE, TN(j+1,i+1)<TE, TLN(j+1,i+1)<TE
        caseID(j,i)=3;
        
        if((TN(j+1,i+1)>TLTM(j+1,i+1))&&(TN(j+1,i+1)>TE)&&(TTM(j+1,i+1)>TLTM(j+1,i+1))&&(FSTTM(j+1,i+1)<=0.0)&&(dQ(j+1,i+1)<0.0)&&(CLCM(j+1,i+1)>=CL0))
            %SITUATION 1: SUPERLIQUIDUS --> SUPERLIQUIDUS COOLING
            caseID(j,i)=1;
        end
        
        if((TN(j+1,i+1)>TLTM(j+1,i+1))&&(TN(j+1,i+1)>TE)&&(TTM(j+1,i+1)>=TLTM(j+1,i+1))&&(FSTTM(j+1,i+1)<=0.0)&&(dQ(j+1,i+1)>0.0)&&(CLCM(j+1,i+1)>=CL0))
            %SITUATION 11: SUPERLIQUIDUS --> SUPERLIQUIDUS HEATING
            caseID(j,i)=11;
        end
        
        if((CLCM(j+1,i+1)>=CE.MgO)&&(CLCM(j+1,i+1)<=CL0)&&(TN(j+1,i+1)>=TE)&&(TN(j+1,i+1)<TLTM(j+1,i+1))&&(FSTTM(j+1,i+1)>=0.0)&&(TTM(j+1,i+1)>=TLTM(j+1,i+1))&&(dQ(j+1,i+1)<0.0))
            %CLCM=MCL.MgO: Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1]
            %SITUATION 2: SUPERLIQUIDUS --> MUSHY REGION COOLING
            caseID(j,i)=2;
        end
        
        if((CLCM(j+1,i+1)>CE.MgO)&&(CLCM(j+1,i+1)<CL0)&&(TN(j+1,i+1)>TE)&&(TN(j+1,i+1)>TLTM(j+1,i+1))&&(TN(j+1,i+1)<=TL0)&&(FSTTM(j+1,i+1)>0.0)&&(TTM(j+1,i+1)>=TE)&&(dQ(j+1,i+1)>0.0)&&(dQ(j+1,i+1)<=QS))
            %CLCM=MCL.MgO: Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1]
            %SITUATION 21: MUSHY REGION --> MUSHY REGION REMELTING
            %local heating but can not heats up to superliquidus
            caseID(j,i)=21;
        end

        if((CLCM(j+1,i+1)>CE.MgO)&&(CLCM(j+1,i+1)<CL0)&&(TN(j+1,i+1)>TE)&&(TN(j+1,i+1)>TLTM(j+1,i+1))&&(TN(j+1,i+1)>TL0)&&(FSTTM(j+1,i+1)>0.0)&&(TTM(j+1,i+1)>=TE)&&(dQ(j+1,i+1)>0.0)&&(dQ(j+1,i+1)>QS))
            %CLCM=MCL.MgO: Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1]
            %SITUATION 22: MUSHY REGION --> SUPERLIQUIDUS REMELTING
            %local heating up to superliquidus
            caseID(j,i)=22;
        end
        
        if((CLCM(j+1,i+1)>CE.MgO)&&(FSTTM(j+1,i+1)>0.0)&&(TN(j+1,i+1)<TE)&&(TTM(j+1,i+1)>TE)&&(dQ(j+1,i+1)<0.0))
            %CLCM=MCL.MgO: Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1]
            %SITUATION 3: MUSHY --> EUTECTIC COOLING
            caseID(j,i)=3;
        end
        
        if((CLCM(j+1,i+1)==CE.MgO)&&(FSTTM(j+1,i+1)>0.0)&&(TN(j+1,i+1)>TE)&&(TTM(j+1,i+1)==TE)&&(dQ(j+1,i+1)>0.0)&&(dQ(j+1,i+1)>=HS.CPX*MS.CPX+HS.PL*MS.PL))
            %CLCM=MCL.MgO: Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1]
            %SITUATION 31: EUTECTIC --> MUSHY HEATING
            %local heating up to mushy region from eutectic
            caseID(j,i)=31;
        end

        if((CLCM(j+1,i+1)<=CE.MgO)&&(FSTTM(j+1,i+1)+dFSLHCAL<=1.0)&&(TN(j+1,i+1)<TE)&&(TTM(j+1,i+1)==TE)&&(dQ(j+1,i+1)<0.0))
            %CLCM=MCL.MgO: Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1]
            %SITUATION 4: EUTECTIC --> EUTECTIC COOLING
            caseID(j,i)=4;
        end
        
        if((CLCM(j+1,i+1)==CE.MgO)&&(FSTTM(j+1,i+1)>0.0)&&(TN(j+1,i+1)>TE)&&(TTM(j+1,i+1)==TE)&&(dQ(j+1,i+1)>0.0)&&(dQ(j+1,i+1)<HS.CPX*MS.CPX+HS.PL*MS.PL))
            %CLCM=MCL.MgO: Critical Major element that represents system evolution, e.g., MgO in liquid of Di-An system [1]
            %SITUATION 41: EUTECTIC --> EUTECTIC HEATING
            %local heating to eutectic from eutectic
            caseID(j,i)=41;
        end
        
        if((CLCM(j+1,i+1)==CE.MgO)&&(TTM(j+1,i+1)==TE)&&(FSTTM(j+1,i+1)<1.0)&&(TN(j+1,i+1)<TE)&&(dFSLHCAL+FSTTM(j+1,i+1)>1.0)&&(dQ(j+1,i+1)<0.0))
            %SITUATION 5: EUTECTIC --> COMPLETE SOLID
            caseID(j,i)=5;
        end
        
        if((CLCM(j+1,i+1)==CE.MgO)&&(TTM(j+1,i+1)<TE)&&(FSTTM(j+1,i+1)==1.0)&&(TN(j+1,i+1)>TE)&&(dQ(j+1,i+1)>0.0)&&(dQ(j+1,i+1)>HS.CPX*MS.CPX+HS.PL*MS.PL))
            %SITUATION 51: COMPLETE SOLID --> EUTECTIC HEATING
            caseID(j,i)=51;
        end
        
        if((CLCM(j+1,i+1)==CE.MgO)&&(TTM(j+1,i+1)<TE)&&(FSTTM(j+1,i+1)==1.0)&&(TN(j+1,i+1)>TE)&&(dQ(j+1,i+1)>0.0)&&(dQ(j+1,i+1)>HS.CPX*MS.CPX+HS.PL*MS.PL+0))
            %SITUATION 52: COMPLETE SOLID --> MUSHY HEATING
            caseID(j,i)=52;
        end
        
        if((FSTTM(j+1,i+1)>=1.0)&&(TTM(j+1,i+1)<=TE)&&(TN(j+1,i+1)<=TE))
            %SITUATION 6: SOLID --> SOLID COOLING & SOLID --> SOLID HEATING
            caseID(j,i)=6;
        end
                
        switch caseID(j,i)
            
            case 1 %SITUATION 1: SUPERLIQUIDUS --> SUPERLIQUIDUS COOLING
                %For superliquidus, we do nothing about dFSLH since no crystals is formed, i.e., dFSLH.ANY=0.0
                dFSLHN.OL(j+1,i+1)=0.0;
                dFSLHN.OPX(j+1,i+1)=0.0;
                dFSLHN.CPX(j+1,i+1)=0.0;
                dFSLHN.PL(j+1,i+1)=0.0;
                dFSLHN.ILM(j+1,i+1)=0.0;
                dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);
                
                CLN.SiO2(j+1,i+1)=MCL.SiO2(j+1,i+1);
                CLN.TiO2(j+1,i+1)=MCL.TiO2(j+1,i+1);
                CLN.Al2O3(j+1,i+1)=MCL.Al2O3(j+1,i+1);
                CLN.FeO(j+1,i+1)=MCL.FeO(j+1,i+1);
                CLN.Fe2O3(j+1,i+1)=MCL.Fe2O3(j+1,i+1);
                CLN.MnO(j+1,i+1)=MCL.MnO(j+1,i+1);
                CLN.MgO(j+1,i+1)=MCL.MgO(j+1,i+1);
                CLN.CaO(j+1,i+1)=MCL.CaO(j+1,i+1);
                CLN.Na2O(j+1,i+1)=MCL.Na2O(j+1,i+1);
                CLN.K2O(j+1,i+1)=MCL.K2O(j+1,i+1);
                CLN.P2O5(j+1,i+1)=MCL.P2O5(j+1,i+1);
                CLN.H2O(j+1,i+1)=MCL.H2O(j+1,i+1);
                CLN.Sm(j+1,i+1)=TCL.Sm(j+1,i+1);
                CLN.Nd(j+1,i+1)=TCL.Nd(j+1,i+1);
                
                iters(j+1,i+1)=1;
                
                LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RS.OL(j+1,i+1)*HS.OL(j+1,i+1);
                LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1);
                LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1);
                LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RS.PL(j+1,i+1)*HS.PL(j+1,i+1);
                LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1);
                LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat [J/m^3]
                
                %latent heat for energy check [J/m^3]
                LAT(j,i)=LHT(j+1,i+1);
                %sensible heat for energy check [J/m^3]
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TN(j+1,i+1);
                %NOTE: RSCPSdFSTM has been calculated in the first guess of TN
                %T for energy check [K]
                TEC(j+1,i+1)=TN(j+1,i+1);
                
                
          case 11 %SITUATION 11: SUPERLIQUIDUS --> SUPERLIQUIDUS HEATING
                %For superliquidus, we do nothing about dFSLH since no crystals is formed, i.e., dFSLH.ANY=0.0
                dFSLHN.OL(j+1,i+1)=0.0;
                dFSLHN.OPX(j+1,i+1)=0.0;
                dFSLHN.CPX(j+1,i+1)=0.0;
                dFSLHN.PL(j+1,i+1)=0.0;
                dFSLHN.ILM(j+1,i+1)=0.0;
                dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);
                
                CLN.SiO2(j+1,i+1)=MCL.SiO2(j+1,i+1);
                CLN.TiO2(j+1,i+1)=MCL.TiO2(j+1,i+1);
                CLN.Al2O3(j+1,i+1)=MCL.Al2O3(j+1,i+1);
                CLN.FeO(j+1,i+1)=MCL.FeO(j+1,i+1);
                CLN.Fe2O3(j+1,i+1)=MCL.Fe2O3(j+1,i+1);
                CLN.MnO(j+1,i+1)=MCL.MnO(j+1,i+1);
                CLN.MgO(j+1,i+1)=MCL.MgO(j+1,i+1);
                CLN.CaO(j+1,i+1)=MCL.CaO(j+1,i+1);
                CLN.Na2O(j+1,i+1)=MCL.Na2O(j+1,i+1);
                CLN.K2O(j+1,i+1)=MCL.K2O(j+1,i+1);
                CLN.P2O5(j+1,i+1)=MCL.P2O5(j+1,i+1);
                CLN.H2O(j+1,i+1)=MCL.H2O(j+1,i+1);
                CLN.Sm(j+1,i+1)=TCL.Sm(j+1,i+1);
                CLN.Nd(j+1,i+1)=TCL.Nd(j+1,i+1);
                
                iters(j+1,i+1)=1;
                
                LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RS.OL(j+1,i+1)*HS.OL(j+1,i+1);
                LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1);
                LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1);
                LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RS.PL(j+1,i+1)*HS.PL(j+1,i+1);
                LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1);
                LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat [J/m^3]
                
                %latent heat for energy check [J/m^3]
                LAT(j,i)=LHT(j+1,i+1);
                %sensible heat for energy check [J/m^3]
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TN(j+1,i+1);
                %NOTE: RSCPSdFSTM has been calculated in the first guess of TN
                %T for energy check [K]
                TEC(j+1,i+1)=TN(j+1,i+1);
                
                
            case 2 %SITUATION 2: SUPERLIQUIDUS --> MUSHY REGION COOLING
                %0.0<=FST<FSE, TE<=TN<=TLTM --> CLCM>CE, TE<=TN<=TLTM
                
                %Once some cell are subliquidus, we mark T has been
                %solved by T-FS-CL iteration. Even one cell!
                modifymarker=modifymarker+1;
                
                %----------------------------- STEP TWO -----------------------------------
                %Take the approximations TN and CLCM(n) as the initial values of temperature and
                %liquid composition for time n+1, and calculate the approximation of CLCM(n+1) using Eq.
                %(9) and the approximate liquidus temperature at n+1 from Eq. (3).
                
                %FL of Next step
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);
                
                %NOTE: new T has effect on RL,RS
                %New RL (=RLN) updated by new T, MCL
                [~,~,RLN(j+1,i+1)]=VisCpRL(TN(j+1,i+1),PATM(j+1,i+1),MCL0);
                
                %New CS updated by new T; originally New KP and then New CS in Xudaming
                [~]=SYSEOS('Di');%MCSN
                
                %New RS (=RSN) updated by new T
                RSNS=SolidDensity(TN(j+1,i+1),PATM(j+1,i+1));%solid phase density [kg/m^3]
                RSN.OL(j+1,i+1)=RSNS.OL;
                RSN.OPX(j+1,i+1)=RSNS.OPX;
                RSN.CPX(j+1,i+1)=RSNS.CPX;
                RSN.PL(j+1,i+1)=RSNS.PL;
                RSN.ILM(j+1,i+1)=RSNS.ILM;
                
                %CLFNS--Factor ahead CL(n+1) for solid part: RSN, MCSN are updated of old RS, MCS; dFSLHTM==dFSLHN==0.0 for the first ste
                %It should be:
                %RSN.Fo(j+1,i+1)*NCSN.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+RS.Fo(j+1,i+1)*NCS.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+...
                %Try TWM, LPUM by pMELTS
                
                CLFNS.SiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.TiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Al2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.FeO(j,i)=RSN.OL(j+1,i+1)*MCSN.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Fe2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.MnO(j,i)=RSN.OL(j+1,i+1)*MCSN.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.MgO(j,i)=RSN.OL(j+1,i+1)*MCSN.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.CaO(j,i)=RSN.OL(j+1,i+1)*MCSN.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Na2O(j,i)=RSN.OL(j+1,i+1)*MCSN.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.K2O(j,i)=RSN.OL(j+1,i+1)*MCSN.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.P2O5(j,i)=RSN.OL(j+1,i+1)*MCSN.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.H2O(j,i)=RSN.OL(j+1,i+1)*MCSN.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                %NOTE: For Trace elements, we use TCSN==TCS as approximation.
                CLFNS.Sm(j,i)=RSN.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Nd(j,i)=RSN.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                %CLFN--Factor ahead CL(n+1) of liquid: FL, RL should use NEW values FLN, RLN
                CLFN.SiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.TiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Al2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.FeO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Fe2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.MnO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.MgO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.CaO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Na2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.K2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.P2O5(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.H2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Sm(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Nd(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                
                
                %New major element CL of Next step in liquid
                CLN.SiO2(j+1,i+1)=(MISS.SiO2(j+1,i+1)+CLQVT.SiO2(j,i)+CLQDT.SiO2(j,i)-CLFNS.SiO2(j,i))/CLFN.SiO2(j,i);%[kg/m^3 divided by kg/m^3 ==1]
                CLN.TiO2(j+1,i+1)=(MISS.TiO2(j+1,i+1)+CLQVT.TiO2(j,i)+CLQDT.TiO2(j,i)-CLFNS.TiO2(j,i))/CLFN.TiO2(j,i);
                CLN.Al2O3(j+1,i+1)=(MISS.Al2O3(j+1,i+1)+CLQVT.Al2O3(j,i)+CLQDT.Al2O3(j,i)-CLFNS.Al2O3(j,i))/CLFN.Al2O3(j,i);
                CLN.FeO(j+1,i+1)=(MISS.FeO(j+1,i+1)+CLQVT.FeO(j,i)+CLQDT.FeO(j,i)-CLFNS.FeO(j,i))/CLFN.FeO(j,i);
                CLN.Fe2O3(j+1,i+1)=(MISS.Fe2O3(j+1,i+1)+CLQVT.Fe2O3(j,i)+CLQDT.Fe2O3(j,i)-CLFNS.Fe2O3(j,i))/CLFN.Fe2O3(j,i);
                CLN.MnO(j+1,i+1)=(MISS.MnO(j+1,i+1)+CLQVT.MnO(j,i)+CLQDT.MnO(j,i)-CLFNS.MnO(j,i))/CLFN.MnO(j,i);
                CLN.MgO(j+1,i+1)=(MISS.MgO(j+1,i+1)+CLQVT.MgO(j,i)+CLQDT.MgO(j,i)-CLFNS.MgO(j,i))/CLFN.MgO(j,i);
                CLN.CaO(j+1,i+1)=(MISS.CaO(j+1,i+1)+CLQVT.CaO(j,i)+CLQDT.CaO(j,i)-CLFNS.CaO(j,i))/CLFN.CaO(j,i);
                CLN.Na2O(j+1,i+1)=(MISS.Na2O(j+1,i+1)+CLQVT.Na2O(j,i)+CLQDT.Na2O(j,i)-CLFNS.Na2O(j,i))/CLFN.Na2O(j,i);
                CLN.K2O(j+1,i+1)=(MISS.K2O(j+1,i+1)+CLQVT.K2O(j,i)+CLQDT.K2O(j,i)-CLFNS.K2O(j,i))/CLFN.K2O(j,i);
                CLN.P2O5(j+1,i+1)=(MISS.P2O5(j+1,i+1)+CLQVT.P2O5(j,i)+CLQDT.P2O5(j,i)-CLFNS.P2O5(j,i))/CLFN.P2O5(j,i);
                CLN.H2O(j+1,i+1)=(MISS.H2O(j+1,i+1)+CLQVT.H2O(j,i)+CLQDT.H2O(j,i)-CLFNS.H2O(j,i))/CLFN.H2O(j,i);
                CLN.Sm(j+1,i+1)=(TISS.Sm(j+1,i+1)+CLQVT.Sm(j,i)+CLQDT.Sm(j,i)-CLFNS.Sm(j,i))/CLFN.Sm(j,i);
                CLN.Nd(j+1,i+1)=(TISS.Nd(j+1,i+1)+CLQVT.Nd(j,i)+CLQDT.Nd(j,i)-CLFNS.Nd(j,i))/CLFN.Nd(j,i);
                
                %Relative error of major elements in liquid [1]
                errCLTM(j+1,i+1)=max([abs((CLN.SiO2(j+1,i+1)-CLTM.SiO2(j+1,i+1))/CLN.SiO2(j+1,i+1)),...
                    abs((CLN.TiO2(j+1,i+1)-CLTM.TiO2(j+1,i+1))/CLN.TiO2(j+1,i+1)),...
                    abs((CLN.Al2O3(j+1,i+1)-CLTM.Al2O3(j+1,i+1))/CLN.Al2O3(j+1,i+1)),...
                    abs((CLN.FeO(j+1,i+1)-CLTM.FeO(j+1,i+1))/CLN.FeO(j+1,i+1)),...
                    abs((CLN.Fe2O3(j+1,i+1)-CLTM.Fe2O3(j+1,i+1))/CLN.Fe2O3(j+1,i+1)),...
                    abs((CLN.MnO(j+1,i+1)-CLTM.MnO(j+1,i+1))/CLN.MnO(j+1,i+1)),...
                    abs((CLN.MgO(j+1,i+1)-CLTM.MgO(j+1,i+1))/CLN.MgO(j+1,i+1)),...
                    abs((CLN.CaO(j+1,i+1)-CLTM.CaO(j+1,i+1))/CLN.CaO(j+1,i+1)),...
                    abs((CLN.Na2O(j+1,i+1)-CLTM.Na2O(j+1,i+1))/CLN.Na2O(j+1,i+1)),...
                    abs((CLN.K2O(j+1,i+1)-CLTM.K2O(j+1,i+1))/CLN.K2O(j+1,i+1)),...
                    abs((CLN.P2O5(j+1,i+1)-CLTM.P2O5(j+1,i+1))/CLN.P2O5(j+1,i+1)),...
                    abs((CLN.H2O(j+1,i+1)-CLTM.H2O(j+1,i+1))/CLN.H2O(j+1,i+1))]);
                
                CLTM.SiO2(j+1,i+1)=CLN.SiO2(j+1,i+1);
                CLTM.TiO2(j+1,i+1)=CLN.TiO2(j+1,i+1);
                CLTM.Al2O3(j+1,i+1)=CLN.Al2O3(j+1,i+1);
                CLTM.FeO(j+1,i+1)=CLN.FeO(j+1,i+1);
                CLTM.Fe2O3(j+1,i+1)=CLN.Fe2O3(j+1,i+1);
                CLTM.MnO(j+1,i+1)=CLN.MnO(j+1,i+1);
                CLTM.MgO(j+1,i+1)=CLN.MgO(j+1,i+1);
                CLTM.CaO(j+1,i+1)=CLN.CaO(j+1,i+1);
                CLTM.Na2O(j+1,i+1)=CLN.Na2O(j+1,i+1);
                CLTM.K2O(j+1,i+1)=CLN.K2O(j+1,i+1);
                CLTM.P2O5(j+1,i+1)=CLN.P2O5(j+1,i+1);
                CLTM.H2O(j+1,i+1)=CLN.H2O(j+1,i+1);
                CLTM.Sm(j+1,i+1)=CLN.Sm(j+1,i+1);
                CLTM.Nd(j+1,i+1)=CLN.Nd(j+1,i+1);
                
                CLCM(j+1,i+1)=CLN.MgO(j+1,i+1);
                
                %Prepare Temporary MCL for VisCpRL.m
                MCL0.SiO2=CLN.SiO2(j+1,i+1);
                MCL0.TiO2=CLN.TiO2(j+1,i+1);
                MCL0.Al2O3=CLN.Al2O3(j+1,i+1);
                MCL0.FeO=CLN.FeO(j+1,i+1);
                MCL0.Fe2O3=CLN.Fe2O3(j+1,i+1);
                MCL0.MnO=CLN.MnO(j+1,i+1);
                MCL0.MgO=CLN.MgO(j+1,i+1);
                MCL0.CaO=CLN.CaO(j+1,i+1);
                MCL0.Na2O=CLN.Na2O(j+1,i+1);
                MCL0.K2O=CLN.K2O(j+1,i+1);
                MCL0.P2O5=CLN.P2O5(j+1,i+1);
                MCL0.H2O=CLN.H2O(j+1,i+1);
                
                %Liquidus is a complex function of liquid composition, TLN=func(MgO, CaO, Al2O3, SiO2, H2O, ...)
                TLN(j+1,i+1)=-1.3351*(CLCM(j+1,i+1)*100.0)^2+55.878*CLCM(j+1,i+1)*100.0+813.9991+273.15;%In Di-An system, new liquidus is a simple function of MgO [K]
                
                iters(j+1,i+1)=0;
                
                m=1;
                %----------------------------- STEP TWO -----------------------------------
                while(err(j+1,i+1))
                    m=m+1;
                    if(m>500)
                        break;
                    end
                    %-------------------------- STEP THREE --------------------------------
                    %Calculate the approximation of dFSLHN by Eq. (10) and take the average:
                    %dFSLHN=0.5*(dFSLHN+dFSLHTM)
                    
                    %STLL: Sensible heat Transform into Latent heat, temporarily used as Left hand side in Eq.(10)
                    %STLR: Sensible heat Transform into Latent heat, temporarily used as Right hand side in Eq.(10)
                    STLL=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLN(j+1,i+1)-TN(j+1,i+1));
                    STLR=RS.OL(j+1,i+1)*HS.OL(j+1,i+1)*0.0+...%dFSLHTM.OL(j+1,i+1)=0.0 in Di-An system
                        RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1)*0.0+...%dFSLHTM.OPX(j+1,i+1)=0.0 in Di-An system
                        RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1)*0.0+...%dFSLHTM.ILM(j+1,i+1)=0.0 in Di-An system
                        RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1)+...
                        RS.PL(j+1,i+1)*HS.PL(j+1,i+1)*0.0;%dFSLHTM.PL(j+1,i+1)=0.0 when PL is not yet liquid phase, i.e., MCL.MgO>CE.MgO
                    
                    dFSLHN.CPX(j+1,i+1)=STLL/STLR;%[1]
                    dFSLHN.OL(j+1,i+1)=0.0;
                    dFSLHN.OPX(j+1,i+1)=0.0;
                    dFSLHN.PL(j+1,i+1)=0.0;
                    dFSLHN.ILM(j+1,i+1)=0.0;
                    
                    dFSLHN.OL(j+1,i+1)=0.5*(dFSLHN.OL(j+1,i+1)+dFSLHTM.OL(j+1,i+1));
                    dFSLHN.OPX(j+1,i+1)=0.5*(dFSLHN.OPX(j+1,i+1)+dFSLHTM.OPX(j+1,i+1));
                    dFSLHN.CPX(j+1,i+1)=0.5*(dFSLHN.CPX(j+1,i+1)+dFSLHTM.CPX(j+1,i+1));
                    dFSLHN.PL(j+1,i+1)=0.5*(dFSLHN.PL(j+1,i+1)+dFSLHTM.PL(j+1,i+1));
                    dFSLHN.ILM(j+1,i+1)=0.5*(dFSLHN.ILM(j+1,i+1)+dFSLHTM.ILM(j+1,i+1));
                    
                    errdFS(j+1,i+1)=max([abs((dFSLHN.OL(j+1,i+1)-dFSLHTM.OL(j+1,i+1))/dFSLHN.OL(j+1,i+1)),...
                        abs((dFSLHN.OPX(j+1,i+1)-dFSLHTM.OPX(j+1,i+1))/dFSLHN.OPX(j+1,i+1)),...
                        abs((dFSLHN.CPX(j+1,i+1)-dFSLHTM.CPX(j+1,i+1))/dFSLHN.CPX(j+1,i+1)),...
                        abs((dFSLHN.PL(j+1,i+1)-dFSLHTM.PL(j+1,i+1))/dFSLHN.PL(j+1,i+1)),...
                        abs((dFSLHN.ILM(j+1,i+1)-dFSLHTM.ILM(j+1,i+1))/dFSLHN.ILM(j+1,i+1))]);
                    
                    dFSLHTM.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1);
                    dFSLHTM.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1);
                    dFSLHTM.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1);
                    dFSLHTM.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1);
                    dFSLHTM.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1);
                    %-------------------------- STEP THREE --------------------------------
                    %                     figure(2)
                    %                     plot(m,TN(j+1,i+1),'g*');
                    %                     hold on
                    %                     plot(m,TLN(j+1,i+1),'r*');
                    %                     %VERY IMPORTANT NOTE: TLN should be >= TN so that dFS >= 0.0, according to
                    %                     %Eq.(10) (sensible heat transformed into latent heat). After this inner
                    %                     %iteration, TN should be TLN bu assuming thermodynamic equilibrium.
                    %                     hold on
                    %                     figure(3)
                    %                     plot(m,CLN.MgO(j+1,i+1),'k*');
                    %                     hold on
                    %                     figure(4)
                    %                     plot(m,dFSLHN.CPX(j+1,i+1),'ko');
                    %                     hold on
                    %                     m=m+1;
                    
                    
                    %--------------------------- STEP FOUR --------------------------------
                    %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
                    %values, dFSLHN and CLN, and take TN=0.5*(TN+TTM).
                    
                    %NOTE: new T, CL have effects on RL, RS, CpS, CpL. For the first guess of T(n+1) and CL(n+1), which are differen from T(n) and CL(n),
                    %should drive updated RL, RS, CpS, CpL immediately. However, according to Temperature Eq., all but dFSLH are old values,
                    %so old RL, RS, CpS, CpS are used, and we don't update them here!!!
                    
                    %latent heat released in action [J/m^3]
                    LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RS.OL(j+1,i+1)*HS.OL(j+1,i+1);
                    LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1);
                    LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1);
                    LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RS.PL(j+1,i+1)*HS.PL(j+1,i+1);
                    LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1);
                    LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat release in action [J/m^3]
                    
                    %Old RS, CpS, New dFSLHN
                    RSCPSdFSTM(j+1,i+1)=RS.OL(j+1,i+1)*CPS.OL(j+1,i+1)*dFSLHN.OL(j+1,i+1)+...
                        RS.OPX(j+1,i+1)*CPS.OPX(j+1,i+1)*dFSLHN.OPX(j+1,i+1)+...
                        RS.CPX(j+1,i+1)*CPS.CPX(j+1,i+1)*dFSLHN.CPX(j+1,i+1)+...
                        RS.PL(j+1,i+1)*CPS.PL(j+1,i+1)*dFSLHN.PL(j+1,i+1)+...
                        RS.ILM(j+1,i+1)*CPS.ILM(j+1,i+1)*dFSLHN.ILM(j+1,i+1);
                    
                    dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                    TFN(j+1,i+1)=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1);%[J/m^3/K]
                    
                    TN(j+1,i+1)=(IHS(j+1,i+1)+TQVT(j,i)+TQDT(j,i)+LHT(j+1,i+1))/TFN(j+1,i+1);%[J/m^3 divided by J/m^3/K == K]
                    TN(j+1,i+1)=0.5*(TN(j+1,i+1)+TTM(j+1,i+1));
                    
                    errTTM(j+1,i+1)=abs((TN(j+1,i+1)-TTM(j+1,i+1))/TN(j+1,i+1));
                    
                    TTM(j+1,i+1)=TN(j+1,i+1);
                    %--------------------------- STEP FOUR --------------------------------
                    
                    %--------------------------- STEP FIVE --------------------------------
                    %With the new initial value TN (as well as dFSLHN and CLN), calculate the new approximation
                    %of CLN(n+1) using Eq. (9), and take CLN=0.5*(CLN+CLTM).
                    
                    %FL of of Next step
                    FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);
                    
                    %NOTE: new T has effect on RL,RS
                    %New RL (=RLN) updated by new T, MCL
                    [~,~,RLN(j+1,i+1)]=VisCpRL(TN(j+1,i+1),PATM(j+1,i+1),MCL0);
                    
                    %New CS updated by new T; originally New KP and then New CS in Xudaming
                    [~]=SYSEOS('Di');%MCSN
                    
                    %New RS (=RSN) updated by new T
                    RSNS=SolidDensity(TN(j+1,i+1),PATM(j+1,i+1));%solid phase density [kg/m^3]
                    RSN.OL(j+1,i+1)=RSNS.OL;
                    RSN.OPX(j+1,i+1)=RSNS.OPX;
                    RSN.CPX(j+1,i+1)=RSNS.CPX;
                    RSN.PL(j+1,i+1)=RSNS.PL;
                    RSN.ILM(j+1,i+1)=RSNS.ILM;
                    
                    %CLFNS--Factor ahead CL(n+1) for solid part: RSN, MCSN are updated of old RS, MCS;
                    %It should be:
                    %RSN.Fo(j+1,i+1)*NCSN.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+RS.Fo(j+1,i+1)*NCS.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+...
                    %Try TWM, LPUM by pMELTS
                    
                    CLFNS.SiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.TiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Al2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.FeO(j,i)=RSN.OL(j+1,i+1)*MCSN.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Fe2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.MnO(j,i)=RSN.OL(j+1,i+1)*MCSN.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.MgO(j,i)=RSN.OL(j+1,i+1)*MCSN.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.CaO(j,i)=RSN.OL(j+1,i+1)*MCSN.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Na2O(j,i)=RSN.OL(j+1,i+1)*MCSN.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.K2O(j,i)=RSN.OL(j+1,i+1)*MCSN.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.P2O5(j,i)=RSN.OL(j+1,i+1)*MCSN.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.H2O(j,i)=RSN.OL(j+1,i+1)*MCSN.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    %NOTE: For Trace elements, we use TCSN==TCS as approximation.
                    CLFNS.Sm(j,i)=RSN.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Nd(j,i)=RSN.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    %CLFN--Factor ahead CL(n+1) of liquid and solid: FL, RL should use NEW values FLN, RLN
                    CLFN.SiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.TiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Al2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.FeO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Fe2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.MnO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.MgO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.CaO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Na2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.K2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.P2O5(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.H2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Sm(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Nd(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    
                    %New major element CL of Next step in liquid [kg/m^3 divided by kg/m^3 ==1]
                    CLN.SiO2(j+1,i+1)=(MISS.SiO2(j+1,i+1)+CLQVT.SiO2(j,i)+CLQDT.SiO2(j,i)-CLFNS.SiO2(j,i))/CLFN.SiO2(j,i);
                    CLN.TiO2(j+1,i+1)=(MISS.TiO2(j+1,i+1)+CLQVT.TiO2(j,i)+CLQDT.TiO2(j,i)-CLFNS.TiO2(j,i))/CLFN.TiO2(j,i);
                    CLN.Al2O3(j+1,i+1)=(MISS.Al2O3(j+1,i+1)+CLQVT.Al2O3(j,i)+CLQDT.Al2O3(j,i)-CLFNS.Al2O3(j,i))/CLFN.Al2O3(j,i);
                    CLN.FeO(j+1,i+1)=(MISS.FeO(j+1,i+1)+CLQVT.FeO(j,i)+CLQDT.FeO(j,i)-CLFNS.FeO(j,i))/CLFN.FeO(j,i);
                    CLN.Fe2O3(j+1,i+1)=(MISS.Fe2O3(j+1,i+1)+CLQVT.Fe2O3(j,i)+CLQDT.Fe2O3(j,i)-CLFNS.Fe2O3(j,i))/CLFN.Fe2O3(j,i);
                    CLN.MnO(j+1,i+1)=(MISS.MnO(j+1,i+1)+CLQVT.MnO(j,i)+CLQDT.MnO(j,i)-CLFNS.MnO(j,i))/CLFN.MnO(j,i);
                    CLN.MgO(j+1,i+1)=(MISS.MgO(j+1,i+1)+CLQVT.MgO(j,i)+CLQDT.MgO(j,i)-CLFNS.MgO(j,i))/CLFN.MgO(j,i);
                    CLN.CaO(j+1,i+1)=(MISS.CaO(j+1,i+1)+CLQVT.CaO(j,i)+CLQDT.CaO(j,i)-CLFNS.CaO(j,i))/CLFN.CaO(j,i);
                    CLN.Na2O(j+1,i+1)=(MISS.Na2O(j+1,i+1)+CLQVT.Na2O(j,i)+CLQDT.Na2O(j,i)-CLFNS.Na2O(j,i))/CLFN.Na2O(j,i);
                    CLN.K2O(j+1,i+1)=(MISS.K2O(j+1,i+1)+CLQVT.K2O(j,i)+CLQDT.K2O(j,i)-CLFNS.K2O(j,i))/CLFN.K2O(j,i);
                    CLN.P2O5(j+1,i+1)=(MISS.P2O5(j+1,i+1)+CLQVT.P2O5(j,i)+CLQDT.P2O5(j,i)-CLFNS.P2O5(j,i))/CLFN.P2O5(j,i);
                    CLN.H2O(j+1,i+1)=(MISS.H2O(j+1,i+1)+CLQVT.H2O(j,i)+CLQDT.H2O(j,i)-CLFNS.H2O(j,i))/CLFN.H2O(j,i);
                    CLN.Sm(j+1,i+1)=(TISS.Sm(j+1,i+1)+CLQVT.Sm(j,i)+CLQDT.Sm(j,i)-CLFNS.Sm(j,i))/CLFN.Sm(j,i);
                    CLN.Nd(j+1,i+1)=(TISS.Nd(j+1,i+1)+CLQVT.Nd(j,i)+CLQDT.Nd(j,i)-CLFNS.Nd(j,i))/CLFN.Nd(j,i);
                    
                    %Take average CLN and CLTM
                    CLN.SiO2(j+1,i+1)=0.5*(CLN.SiO2(j+1,i+1)+CLTM.SiO2(j+1,i+1));
                    CLN.TiO2(j+1,i+1)=0.5*(CLN.TiO2(j+1,i+1)+CLTM.TiO2(j+1,i+1));
                    CLN.Al2O3(j+1,i+1)=0.5*(CLN.Al2O3(j+1,i+1)+CLTM.Al2O3(j+1,i+1));
                    CLN.FeO(j+1,i+1)=0.5*(CLN.FeO(j+1,i+1)+CLTM.FeO(j+1,i+1));
                    CLN.Fe2O3(j+1,i+1)=0.5*(CLN.Fe2O3(j+1,i+1)+CLTM.Fe2O3(j+1,i+1));
                    CLN.MnO(j+1,i+1)=0.5*(CLN.MnO(j+1,i+1)+CLTM.MnO(j+1,i+1));
                    CLN.MgO(j+1,i+1)=0.5*(CLN.MgO(j+1,i+1)+CLTM.MgO(j+1,i+1));
                    CLN.CaO(j+1,i+1)=0.5*(CLN.CaO(j+1,i+1)+CLTM.CaO(j+1,i+1));
                    CLN.Na2O(j+1,i+1)=0.5*(CLN.Na2O(j+1,i+1)+CLTM.Na2O(j+1,i+1));
                    CLN.K2O(j+1,i+1)=0.5*(CLN.K2O(j+1,i+1)+CLTM.K2O(j+1,i+1));
                    CLN.P2O5(j+1,i+1)=0.5*(CLN.P2O5(j+1,i+1)+CLTM.P2O5(j+1,i+1));
                    CLN.H2O(j+1,i+1)=0.5*(CLN.H2O(j+1,i+1)+CLTM.H2O(j+1,i+1));
                    CLN.Sm(j+1,i+1)=0.5*(CLN.Sm(j+1,i+1)+CLTM.Sm(j+1,i+1));
                    CLN.Nd(j+1,i+1)=0.5*(CLN.Nd(j+1,i+1)+CLTM.Nd(j+1,i+1));
                    
                    %Relative error of major elements in liquid [1]
                    errCLTM(j+1,i+1)=max([abs((CLN.SiO2(j+1,i+1)-CLTM.SiO2(j+1,i+1))/CLN.SiO2(j+1,i+1)),...
                        abs((CLN.TiO2(j+1,i+1)-CLTM.TiO2(j+1,i+1))/CLN.TiO2(j+1,i+1)),...
                        abs((CLN.Al2O3(j+1,i+1)-CLTM.Al2O3(j+1,i+1))/CLN.Al2O3(j+1,i+1)),...
                        abs((CLN.FeO(j+1,i+1)-CLTM.FeO(j+1,i+1))/CLN.FeO(j+1,i+1)),...
                        abs((CLN.Fe2O3(j+1,i+1)-CLTM.Fe2O3(j+1,i+1))/CLN.Fe2O3(j+1,i+1)),...
                        abs((CLN.MnO(j+1,i+1)-CLTM.MnO(j+1,i+1))/CLN.MnO(j+1,i+1)),...
                        abs((CLN.MgO(j+1,i+1)-CLTM.MgO(j+1,i+1))/CLN.MgO(j+1,i+1)),...
                        abs((CLN.CaO(j+1,i+1)-CLTM.CaO(j+1,i+1))/CLN.CaO(j+1,i+1)),...
                        abs((CLN.Na2O(j+1,i+1)-CLTM.Na2O(j+1,i+1))/CLN.Na2O(j+1,i+1)),...
                        abs((CLN.K2O(j+1,i+1)-CLTM.K2O(j+1,i+1))/CLN.K2O(j+1,i+1)),...
                        abs((CLN.P2O5(j+1,i+1)-CLTM.P2O5(j+1,i+1))/CLN.P2O5(j+1,i+1)),...
                        abs((CLN.H2O(j+1,i+1)-CLTM.H2O(j+1,i+1))/CLN.H2O(j+1,i+1))]);
                    
                    %Exchange old and new CL
                    CLTM.SiO2(j+1,i+1)=CLN.SiO2(j+1,i+1);
                    CLTM.TiO2(j+1,i+1)=CLN.TiO2(j+1,i+1);
                    CLTM.Al2O3(j+1,i+1)=CLN.Al2O3(j+1,i+1);
                    CLTM.FeO(j+1,i+1)=CLN.FeO(j+1,i+1);
                    CLTM.Fe2O3(j+1,i+1)=CLN.Fe2O3(j+1,i+1);
                    CLTM.MnO(j+1,i+1)=CLN.MnO(j+1,i+1);
                    CLTM.MgO(j+1,i+1)=CLN.MgO(j+1,i+1);
                    CLTM.CaO(j+1,i+1)=CLN.CaO(j+1,i+1);
                    CLTM.Na2O(j+1,i+1)=CLN.Na2O(j+1,i+1);
                    CLTM.K2O(j+1,i+1)=CLN.K2O(j+1,i+1);
                    CLTM.P2O5(j+1,i+1)=CLN.P2O5(j+1,i+1);
                    CLTM.H2O(j+1,i+1)=CLN.H2O(j+1,i+1);
                    CLTM.Sm(j+1,i+1)=CLN.Sm(j+1,i+1);
                    CLTM.Nd(j+1,i+1)=CLN.Nd(j+1,i+1);
                    
                    %--------------------------- STEP FIVE --------------------------------
                    
                    %Prepare Temporary MCL for VisCpRL.m
                    MCL0.SiO2=CLN.SiO2(j+1,i+1);
                    MCL0.TiO2=CLN.TiO2(j+1,i+1);
                    MCL0.Al2O3=CLN.Al2O3(j+1,i+1);
                    MCL0.FeO=CLN.FeO(j+1,i+1);
                    MCL0.Fe2O3=CLN.Fe2O3(j+1,i+1);
                    MCL0.MnO=CLN.MnO(j+1,i+1);
                    MCL0.MgO=CLN.MgO(j+1,i+1);
                    MCL0.CaO=CLN.CaO(j+1,i+1);
                    MCL0.Na2O=CLN.Na2O(j+1,i+1);
                    MCL0.K2O=CLN.K2O(j+1,i+1);
                    MCL0.P2O5=CLN.P2O5(j+1,i+1);
                    MCL0.H2O=CLN.H2O(j+1,i+1);
                    
                    %--------------------------- STEP SIX ---------------------------------
                    err(j+1,i+1)=logical((err0dFSLH<=errdFS(j+1,i+1))||(err0TTM<=errTTM(j+1,i+1))||(err0CLTM<=errCLTM(j+1,i+1)));%logical variable for while loop
                    
                    %If the conditions are not true, take TLN(n+1)=TLN(CLN(n+1)) and repeat steps 3-6 until
                    %all the above conditions are satisfied
                    
                    %--------------------------- STEP SIX ---------------------------------
                    
                    %--------------------------- STEP SEVEN -------------------------------
                    %If jugments are not true, then repeat step 3-6 with updated TN, CLN, dFSLHN
                    
                    %Critical major element taht represents system evolution
                    CLCM(j+1,i+1)=CLN.MgO(j+1,i+1);
                    
                    %Liquidus is a complex function of liquid composition, TLN=func(MgO, CaO, Al2O3, SiO2, H2O, ...)
                    TLN(j+1,i+1)=-1.3351*(CLCM(j+1,i+1)*100.0)^2+55.878*CLCM(j+1,i+1)*100.0+813.9991+273.15;%In Di-An system, new liquidus is a simple function of MgO [K]
                    
                    %-------------------------- STEP SEVEN --------------------------------
                    iters(j+1,i+1)=iters(j+1,i+1)+1;
                end
                
                %latent heat for energy check [J/m^3]
                LAT(j,i)=LHT(j+1,i+1);
                %sensible heat for energy check [J/m^3]
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TN(j+1,i+1);
                %T for energy check [K]
                TEC(j+1,i+1)=TN(j+1,i+1);
                
                %THERMAL EQUILIBRIUM --> T==Tliq
                TN(j+1,i+1)=TLN(j+1,i+1);
                
                
          case 21 %SITUATION 11: MUSHY --> MUSHY HEATING
                
                dFSLHN.OL(j+1,i+1)=0.0;
                dFSLHN.OPX(j+1,i+1)=0.0;
                dFSLHN.CPX(j+1,i+1)=0;
                dFSLHN.PL(j+1,i+1)=0.0;
                dFSLHN.ILM(j+1,i+1)=0.0;
                dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);
                
                CLN.SiO2(j+1,i+1)=MCL.SiO2(j+1,i+1);
                CLN.TiO2(j+1,i+1)=MCL.TiO2(j+1,i+1);
                CLN.Al2O3(j+1,i+1)=MCL.Al2O3(j+1,i+1);
                CLN.FeO(j+1,i+1)=MCL.FeO(j+1,i+1);
                CLN.Fe2O3(j+1,i+1)=MCL.Fe2O3(j+1,i+1);
                CLN.MnO(j+1,i+1)=MCL.MnO(j+1,i+1);
                CLN.MgO(j+1,i+1)=MCL.MgO(j+1,i+1);
                CLN.CaO(j+1,i+1)=MCL.CaO(j+1,i+1);
                CLN.Na2O(j+1,i+1)=MCL.Na2O(j+1,i+1);
                CLN.K2O(j+1,i+1)=MCL.K2O(j+1,i+1);
                CLN.P2O5(j+1,i+1)=MCL.P2O5(j+1,i+1);
                CLN.H2O(j+1,i+1)=MCL.H2O(j+1,i+1);
                CLN.Sm(j+1,i+1)=TCL.Sm(j+1,i+1);
                CLN.Nd(j+1,i+1)=TCL.Nd(j+1,i+1);
                
                iters(j+1,i+1)=1;
                
                LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RS.OL(j+1,i+1)*HS.OL(j+1,i+1);
                LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1);
                LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1);
                LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RS.PL(j+1,i+1)*HS.PL(j+1,i+1);
                LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1);
                LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat [J/m^3]
                
                %latent heat for energy check [J/m^3]
                LAT(j,i)=LHT(j+1,i+1);
                %sensible heat for energy check [J/m^3]
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TN(j+1,i+1);
                %NOTE: RSCPSdFSTM has been calculated in the first guess of TN
                %T for energy check [K]
                TEC(j+1,i+1)=TN(j+1,i+1);
                
                
            case 4 %SITUATION 4: EUTECTIC
                modifymarker=modifymarker+1;
                
                dFSLHN.OL(j+1,i+1)=0.0;
                dFSLHN.OPX(j+1,i+1)=0.0;
                dFSLHN.CPX(j+1,i+1)=dFSLHCPX;
                dFSLHN.PL(j+1,i+1)=dFSLHPL;
                dFSLHN.ILM(j+1,i+1)=0.0;
                dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                
                LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RSE.OL*HS.OL(j+1,i+1);
                LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RSE.OPX*HS.OPX(j+1,i+1);
                LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RSE.CPX*HS.CPX(j+1,i+1);
                LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RSE.PL*HS.PL(j+1,i+1);
                LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RSE.ILM*HS.ILM(j+1,i+1);
                LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat [J/m^3]
                
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);%dFSSMT(j,i) ->0.0 since FSTTM ->1.0
                
                TN(j+1,i+1)=TE;
                
                CLN.SiO2(j+1,i+1)=CE.SiO2;
                CLN.TiO2(j+1,i+1)=CE.TiO2;
                CLN.Al2O3(j+1,i+1)=CE.Al2O3;
                CLN.FeO(j+1,i+1)=CE.FeO;
                CLN.Fe2O3(j+1,i+1)=CE.Fe2O3;
                CLN.MnO(j+1,i+1)=CE.MnO;
                CLN.MgO(j+1,i+1)=CE.MgO;
                CLN.CaO(j+1,i+1)=CE.CaO;
                CLN.Na2O(j+1,i+1)=CE.Na2O;
                CLN.K2O(j+1,i+1)=CE.K2O;
                CLN.P2O5(j+1,i+1)=CE.P2O5;
                CLN.H2O(j+1,i+1)=CE.H2O;
                
                %NOTE: We use last step TCL as eutectic TCL for simplicity.
                CLN.Sm(j+1,i+1)=TCL.Sm(j+1,i+1);
                CLN.Nd(j+1,i+1)=TCL.Nd(j+1,i+1);
                
                iters(j+1,i+1)=1;
                
                TEC(j+1,i+1)=TN(j+1,i+1);%T for energy check
                
                %latent heat for energy check [J/m^3]
                LAT(j,i)=LHT(j+1,i+1);
                
                %sensible heat for energy check [J/m^3]
                RSCPSdFSTM(j+1,i+1)=RSE.OL*CPSE.OL*dFSLHN.OL(j+1,i+1)+...
                    RSE.OPX*CPSE.OPX*dFSLHN.OPX(j+1,i+1)+...
                    RSE.CPX*CPSE.CPX*dFSLHN.CPX(j+1,i+1)+...
                    RSE.PL*CPSE.PL*dFSLHN.PL(j+1,i+1)+...
                    RSE.ILM*CPSE.ILM*dFSLHN.ILM(j+1,i+1);%Old RS, CpS, New dFSLHN
                
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TN(j+1,i+1);
                %NOTE: RLTM==RLE, CPLTM==CPLE, TN==TE==TTMBE
                
                %NOTE: dFSLHN.ANY(j+1,i+1)*(RSCPSTM.ANY(j+1,i+1)-RLCPLTM(j+1,i+1))*TE is sensible heat
                
            case 5 %SITUATION 5: EUTECTIC --> COMPLETE SOLID
                modifymarker=modifymarker+1;
                
                dFSLHT(j+1,i+1)=1.0-FSTTM(j+1,i+1);
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);%here, dFSSMT(j,i)=0.0 since FSTTM ->1.0
                
                dFSLHN.OL(j+1,i+1)=0.0;
                dFSLHN.OPX(j+1,i+1)=0.0;
                dFSLHN.PL(j+1,i+1)=dFSLHT(j+1,i+1)/(1.0+0.599*RSE.PL/(0.401*RSE.CPX));
                dFSLHN.CPX(j+1,i+1)=dFSLHN.PL(j+1,i+1)*0.599*RSE.PL/(0.401*RSE.CPX);
                dFSLHN.ILM(j+1,i+1)=0.0;
                
                % TN(j+1,i+1)=TN(j+1,i+1)+(0.599*RSE.PL*HS.CPX(j+1,i+1)/0.401+RSE.PL*HS.PL(j+1,i+1)+...%latent heat part
                % 0.599*RSE.PL*CPSE.CPX*TE/0.401+RSE.PL*CPSE.PL*TE-(1.0+0.599*RSE.PL/(0.401*RSE.CPX))*RLTM(j+1,i+1)*CPLTM(j+1,i+1)*TE)*(1.0-FSTTM(j+1,i+1))/...%solid property change sensible heat
                % (FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1));
                % TN(j+1,i+1)=TN(j+1,i+1)+(dFSLHN.PL(j+1,i+1)*RSE.PL*HS.PL(j+1,i+1)+dFSLHN.CPX(j+1,i+1)*RSE.CPX*HS.CPX(j+1,i+1)...%latent heat
                %     -dFSLHN.PL(j+1,i+1)*RSE.PL*CPSE.PL*TE-dFSLHN.CPX(j+1,i+1)*RSE.CPX*CPSE.CPX*TE+RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1)*TE)...%solid property change
                % /(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1));%RLTM==RLE, CPLTM==CPLE
                
                %TN(j+1,i+1)=TN(j+1,i+1)+(RSE*HSE-(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*TE)*(1.0-FSTM(j+1,i+1))/(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1));%Eq.(15) in Xudaming
                
                CLN.SiO2(j+1,i+1)=CE.SiO2;
                CLN.TiO2(j+1,i+1)=CE.TiO2;
                CLN.Al2O3(j+1,i+1)=CE.Al2O3;
                CLN.FeO(j+1,i+1)=CE.FeO;
                CLN.Fe2O3(j+1,i+1)=CE.Fe2O3;
                CLN.MnO(j+1,i+1)=CE.MnO;
                CLN.MgO(j+1,i+1)=CE.MgO;
                CLN.CaO(j+1,i+1)=CE.CaO;
                CLN.Na2O(j+1,i+1)=CE.Na2O;
                CLN.K2O(j+1,i+1)=CE.K2O;
                CLN.P2O5(j+1,i+1)=CE.P2O5;
                CLN.H2O(j+1,i+1)=CE.H2O;
                %NOTE: We use last step TCL as eutectic TCL for simplicity.
                CLN.Sm(j+1,i+1)=TCL.Sm(j+1,i+1);
                CLN.Nd(j+1,i+1)=TCL.Nd(j+1,i+1);
                
                LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RSE.OL*HS.OL(j+1,i+1);
                LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RSE.OPX*HS.OPX(j+1,i+1);
                LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RSE.CPX*HS.CPX(j+1,i+1);
                LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RSE.PL*HS.PL(j+1,i+1);
                LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RSE.ILM*HS.ILM(j+1,i+1);
                LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat [J/m^3]
                
                TN(j+1,i+1)=TE-((FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TTMBE(j+1,i+1)-TN(j+1,i+1))-LHT(j+1,i+1)+...
                    dFSLHN.PL(j+1,i+1)*RSE.PL*CPSE.PL*TE+dFSLHN.CPX(j+1,i+1)*RSE.CPX*CPSE.CPX*TE-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1)*TE)/...
                    (FRCPSTM(j+1,i+1)+dFSLHN.PL(j+1,i+1)*RSE.PL*CPSE.PL+dFSLHN.CPX(j+1,i+1)*RSE.CPX*CPSE.CPX);
                
                %latent heat for energy check [J/m^3]
                LAT(j,i)=LHT(j+1,i+1);
                
                %sensible heat for energy check [J/m^3]
                RSCPSdFSTM(j+1,i+1)=RSE.OL*CPSE.OL*dFSLHN.OL(j+1,i+1)+...
                    RSE.OPX*CPSE.OPX*dFSLHN.OPX(j+1,i+1)+...
                    RSE.CPX*CPSE.CPX*dFSLHN.CPX(j+1,i+1)+...
                    RSE.PL*CPSE.PL*dFSLHN.PL(j+1,i+1)+...
                    RSE.ILM*CPSE.ILM*dFSLHN.ILM(j+1,i+1);%Old RS, CpS, New dFSLHN
                
                SEN(j,i)=(FRCPSTM(j+1,i+1)+dFSLHN.PL(j+1,i+1)*RSE.PL*CPSE.PL+dFSLHN.CPX(j+1,i+1)*RSE.CPX*CPSE.CPX)*(TN(j+1,i+1)-TE)+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TE;
                %NOTE: RLTM==RLE, CPLTM==CPLE, TTMBE==TE
                
                %T for energy check [K]
                TEC(j+1,i+1)=TN(j+1,i+1);
                
                iters(j+1,i+1)=1;
                
            case 6 %COMPLETE SOLID COOLING
                
                FLN(j+1,i+1)=0.0;
                
                dFSLHN.OL(j+1,i+1)=0.0;
                dFSLHN.OPX(j+1,i+1)=0.0;
                dFSLHN.PL(j+1,i+1)=0.0;
                dFSLHN.CPX(j+1,i+1)=0.0;
                dFSLHN.ILM(j+1,i+1)=0.0;
                dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                
                CLN.SiO2(j+1,i+1)=CE.SiO2;
                CLN.TiO2(j+1,i+1)=CE.TiO2;
                CLN.Al2O3(j+1,i+1)=CE.Al2O3;
                CLN.FeO(j+1,i+1)=CE.FeO;
                CLN.Fe2O3(j+1,i+1)=CE.Fe2O3;
                CLN.MnO(j+1,i+1)=CE.MnO;
                CLN.MgO(j+1,i+1)=CE.MgO;
                CLN.CaO(j+1,i+1)=CE.CaO;
                CLN.Na2O(j+1,i+1)=CE.Na2O;
                CLN.K2O(j+1,i+1)=CE.K2O;
                CLN.P2O5(j+1,i+1)=CE.P2O5;
                CLN.H2O(j+1,i+1)=CE.H2O;
                %NOTE: We use last step TCL as eutectic TCL for simplicity.
                CLN.Sm(j+1,i+1)=TCL.Sm(j+1,i+1);
                CLN.Nd(j+1,i+1)=TCL.Nd(j+1,i+1);
                
                %First guess of TN(n+1) is just right
                %TN(j+1,i+1)=...
                
                LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RSE.OL*HS.OL(j+1,i+1);
                LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RSE.OPX*HS.OPX(j+1,i+1);
                LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RSE.CPX*HS.CPX(j+1,i+1);
                LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RSE.PL*HS.PL(j+1,i+1);
                LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RSE.ILM*HS.ILM(j+1,i+1);
                LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat [J/m^3]
                
                %latent heat for energy check [J/m^3]
                LAT(j,i)=LHT(j+1,i+1);
                
                %sensible heat for energy check [J/m^3]
                RSCPSdFSTM(j+1,i+1)=RSE.OL*CPSE.OL*dFSLHN.OL(j+1,i+1)+...
                    RSE.OPX*CPSE.OPX*dFSLHN.OPX(j+1,i+1)+...
                    RSE.CPX*CPSE.CPX*dFSLHN.CPX(j+1,i+1)+...
                    RSE.PL*CPSE.PL*dFSLHN.PL(j+1,i+1)+...
                    RSE.ILM*CPSE.ILM*dFSLHN.ILM(j+1,i+1);%Old RS, CpS, New dFSLHN
                
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TN(j+1,i+1);
                %NOTE: RLTM==RLE, CPLTM==CPLE, FLTM=0.0, TTMBE<TE, RSCPSdFSTM==0.0, dFSLHT==0.0
                
                %T for energy check [K]
                TEC(j+1,i+1)=TN(j+1,i+1);
                
                iters(j+1,i+1)=1;
                
            otherwise %SITUATION 3: MUSHY REGION --> EUTECTIC
                
                %                 TNCAL=TN(j+1,i+1);
                %                 dFSTMCAL=dFSTM(j+1,i+1);
                %                 dFSLHTMCAL.OL=dFSLHTM.OL(j+1,i+1);
                %
                %                 TTMCAL=TTM(j+1,i+1);
                %                 %FL of Next step
                %                 FLNCAL=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);%dFSSMT(j,i) probably ->0.0 since it's going to eutectic point
                %If(FLNCAL>0.0) correspoding to statement in the following: fprintf('Approaching eutectic\n');
                
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);
                %NOTE: new T has effect on RL, RS
                %New RL (=RLN) updated by new T, MCL
                [~,~,RLN(j+1,i+1)]=VisCpRL(TN(j+1,i+1),PATM(j+1,i+1),MCL0);
                
                %New CS updated by new T; originally New KP and then New CS in Xudaming
                [~]=SYSEOS('Di');%MCSN
                
                %New RS (=RSN) updated by new T
                RSNS=SolidDensity(TN(j+1,i+1),PATM(j+1,i+1));%solid phase density [kg/m^3]
                RSN.OL(j+1,i+1)=RSNS.OL;
                RSN.OPX(j+1,i+1)=RSNS.OPX;
                RSN.CPX(j+1,i+1)=RSNS.CPX;
                RSN.PL(j+1,i+1)=RSNS.PL;
                RSN.ILM(j+1,i+1)=RSNS.ILM;
                
                %CLFNS--Factor ahead CL(n+1) for solid part: RSN, MCSN are updated of old RS, MCS; dFSLHTM==dFSLHN==0.0 for the first ste
                %It should be:
                %RSN.Fo(j+1,i+1)*NCSN.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+RS.Fo(j+1,i+1)*NCS.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+...
                %Try TWM, LPUM by pMELTS
                
                CLFNS.SiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.TiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Al2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.FeO(j,i)=RSN.OL(j+1,i+1)*MCSN.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Fe2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.MnO(j,i)=RSN.OL(j+1,i+1)*MCSN.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.MgO(j,i)=RSN.OL(j+1,i+1)*MCSN.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.CaO(j,i)=RSN.OL(j+1,i+1)*MCSN.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Na2O(j,i)=RSN.OL(j+1,i+1)*MCSN.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.K2O(j,i)=RSN.OL(j+1,i+1)*MCSN.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.P2O5(j,i)=RSN.OL(j+1,i+1)*MCSN.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.H2O(j,i)=RSN.OL(j+1,i+1)*MCSN.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*MCSN.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*MCSN.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*MCSN.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*MCSN.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                %NOTE: For Trace elements, we use TCSN==TCS as approximation.
                CLFNS.Sm(j,i)=RSN.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                CLFNS.Nd(j,i)=RSN.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                    RSN.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                    RSN.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                    RSN.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                    RSN.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                
                %CLFN--Factor ahead CL(n+1) of liquid: FL, RL should use NEW values FLN, RLN
                CLFN.SiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.TiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Al2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.FeO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Fe2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.MnO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.MgO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.CaO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Na2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.K2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.P2O5(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.H2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Sm(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                CLFN.Nd(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                
                %New major element CL of Next step in liquid
                CLN.SiO2(j+1,i+1)=(MISS.SiO2(j+1,i+1)+CLQVT.SiO2(j,i)+CLQDT.SiO2(j,i)-CLFNS.SiO2(j,i))/CLFN.SiO2(j,i);%[kg/m^3 divided by kg/m^3 ==1]
                CLN.TiO2(j+1,i+1)=(MISS.TiO2(j+1,i+1)+CLQVT.TiO2(j,i)+CLQDT.TiO2(j,i)-CLFNS.TiO2(j,i))/CLFN.TiO2(j,i);
                CLN.Al2O3(j+1,i+1)=(MISS.Al2O3(j+1,i+1)+CLQVT.Al2O3(j,i)+CLQDT.Al2O3(j,i)-CLFNS.Al2O3(j,i))/CLFN.Al2O3(j,i);
                CLN.FeO(j+1,i+1)=(MISS.FeO(j+1,i+1)+CLQVT.FeO(j,i)+CLQDT.FeO(j,i)-CLFNS.FeO(j,i))/CLFN.FeO(j,i);
                CLN.Fe2O3(j+1,i+1)=(MISS.Fe2O3(j+1,i+1)+CLQVT.Fe2O3(j,i)+CLQDT.Fe2O3(j,i)-CLFNS.Fe2O3(j,i))/CLFN.Fe2O3(j,i);
                CLN.MnO(j+1,i+1)=(MISS.MnO(j+1,i+1)+CLQVT.MnO(j,i)+CLQDT.MnO(j,i)-CLFNS.MnO(j,i))/CLFN.MnO(j,i);
                CLN.MgO(j+1,i+1)=(MISS.MgO(j+1,i+1)+CLQVT.MgO(j,i)+CLQDT.MgO(j,i)-CLFNS.MgO(j,i))/CLFN.MgO(j,i);
                CLN.CaO(j+1,i+1)=(MISS.CaO(j+1,i+1)+CLQVT.CaO(j,i)+CLQDT.CaO(j,i)-CLFNS.CaO(j,i))/CLFN.CaO(j,i);
                CLN.Na2O(j+1,i+1)=(MISS.Na2O(j+1,i+1)+CLQVT.Na2O(j,i)+CLQDT.Na2O(j,i)-CLFNS.Na2O(j,i))/CLFN.Na2O(j,i);
                CLN.K2O(j+1,i+1)=(MISS.K2O(j+1,i+1)+CLQVT.K2O(j,i)+CLQDT.K2O(j,i)-CLFNS.K2O(j,i))/CLFN.K2O(j,i);
                CLN.P2O5(j+1,i+1)=(MISS.P2O5(j+1,i+1)+CLQVT.P2O5(j,i)+CLQDT.P2O5(j,i)-CLFNS.P2O5(j,i))/CLFN.P2O5(j,i);
                CLN.H2O(j+1,i+1)=(MISS.H2O(j+1,i+1)+CLQVT.H2O(j,i)+CLQDT.H2O(j,i)-CLFNS.H2O(j,i))/CLFN.H2O(j,i);
                CLN.Sm(j+1,i+1)=(TISS.Sm(j+1,i+1)+CLQVT.Sm(j,i)+CLQDT.Sm(j,i)-CLFNS.Sm(j,i))/CLFN.Sm(j,i);
                CLN.Nd(j+1,i+1)=(TISS.Nd(j+1,i+1)+CLQVT.Nd(j,i)+CLQDT.Nd(j,i)-CLFNS.Nd(j,i))/CLFN.Nd(j,i);
                
                %Relative error of major elements in liquid [1]
                errCLTM(j+1,i+1)=max([abs((CLN.SiO2(j+1,i+1)-CLTM.SiO2(j+1,i+1))/CLN.SiO2(j+1,i+1)),...
                    abs((CLN.TiO2(j+1,i+1)-CLTM.TiO2(j+1,i+1))/CLN.TiO2(j+1,i+1)),...
                    abs((CLN.Al2O3(j+1,i+1)-CLTM.Al2O3(j+1,i+1))/CLN.Al2O3(j+1,i+1)),...
                    abs((CLN.FeO(j+1,i+1)-CLTM.FeO(j+1,i+1))/CLN.FeO(j+1,i+1)),...
                    abs((CLN.Fe2O3(j+1,i+1)-CLTM.Fe2O3(j+1,i+1))/CLN.Fe2O3(j+1,i+1)),...
                    abs((CLN.MnO(j+1,i+1)-CLTM.MnO(j+1,i+1))/CLN.MnO(j+1,i+1)),...
                    abs((CLN.MgO(j+1,i+1)-CLTM.MgO(j+1,i+1))/CLN.MgO(j+1,i+1)),...
                    abs((CLN.CaO(j+1,i+1)-CLTM.CaO(j+1,i+1))/CLN.CaO(j+1,i+1)),...
                    abs((CLN.Na2O(j+1,i+1)-CLTM.Na2O(j+1,i+1))/CLN.Na2O(j+1,i+1)),...
                    abs((CLN.K2O(j+1,i+1)-CLTM.K2O(j+1,i+1))/CLN.K2O(j+1,i+1)),...
                    abs((CLN.P2O5(j+1,i+1)-CLTM.P2O5(j+1,i+1))/CLN.P2O5(j+1,i+1)),...
                    abs((CLN.H2O(j+1,i+1)-CLTM.H2O(j+1,i+1))/CLN.H2O(j+1,i+1))]);
                
                CLTM.SiO2(j+1,i+1)=CLN.SiO2(j+1,i+1);
                CLTM.TiO2(j+1,i+1)=CLN.TiO2(j+1,i+1);
                CLTM.Al2O3(j+1,i+1)=CLN.Al2O3(j+1,i+1);
                CLTM.FeO(j+1,i+1)=CLN.FeO(j+1,i+1);
                CLTM.Fe2O3(j+1,i+1)=CLN.Fe2O3(j+1,i+1);
                CLTM.MnO(j+1,i+1)=CLN.MnO(j+1,i+1);
                CLTM.MgO(j+1,i+1)=CLN.MgO(j+1,i+1);
                CLTM.CaO(j+1,i+1)=CLN.CaO(j+1,i+1);
                CLTM.Na2O(j+1,i+1)=CLN.Na2O(j+1,i+1);
                CLTM.K2O(j+1,i+1)=CLN.K2O(j+1,i+1);
                CLTM.P2O5(j+1,i+1)=CLN.P2O5(j+1,i+1);
                CLTM.H2O(j+1,i+1)=CLN.H2O(j+1,i+1);
                CLTM.Sm(j+1,i+1)=CLN.Sm(j+1,i+1);
                CLTM.Nd(j+1,i+1)=CLN.Nd(j+1,i+1);
                
                %Prepare Temporary MCL for VisCpRL.m
                MCL0.SiO2=CLN.SiO2(j+1,i+1);
                MCL0.TiO2=CLN.TiO2(j+1,i+1);
                MCL0.Al2O3=CLN.Al2O3(j+1,i+1);
                MCL0.FeO=CLN.FeO(j+1,i+1);
                MCL0.Fe2O3=CLN.Fe2O3(j+1,i+1);
                MCL0.MnO=CLN.MnO(j+1,i+1);
                MCL0.MgO=CLN.MgO(j+1,i+1);
                MCL0.CaO=CLN.CaO(j+1,i+1);
                MCL0.Na2O=CLN.Na2O(j+1,i+1);
                MCL0.K2O=CLN.K2O(j+1,i+1);
                MCL0.P2O5=CLN.P2O5(j+1,i+1);
                MCL0.H2O=CLN.H2O(j+1,i+1);
                
                CLCM(j+1,i+1)=CLN.MgO(j+1,i+1);
                
                %Liquidus is a complex function of liquid composition, TLN=func(MgO, CaO, Al2O3, SiO2, H2O, ...)
                TLN(j+1,i+1)=-1.3351*(CLCM(j+1,i+1)*100.0)^2+55.878*CLCM(j+1,i+1)*100.0+813.9991+273.15;%In Di-An system, new liquidus is a simple function of MgO [K]
                
                iters(j+1,i+1)=0;
                m=1;
                %----------------------------- STEP TWO -----------------------------------
                while(err(j+1,i+1))
                    m=m+1;
                    if(m>500)
                        break;
                    end
                    %-------------------------- STEP THREE --------------------------------
                    %Calculate the approximation of dFSLHN by Eq. (10) and take the average:
                    %dFSLHN=0.5*(dFSLHN+dFSLHTM)
                    
                    %STLL: Sensible heat Transform into Latent heat, temporarily used as Left hand side in Eq.(10)
                    %STLR: Sensible heat Transform into Latent heat, temporarily used as Right hand side in Eq.(10)
                    STLL=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLN(j+1,i+1)-TN(j+1,i+1));
                    STLR=RS.OL(j+1,i+1)*HS.OL(j+1,i+1)*0.0+...%dFSLHTM.OL(j+1,i+1)=0.0 in Di-An system
                        RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1)*0.0+...%dFSLHTM.OPX(j+1,i+1)=0.0 in Di-An system
                        RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1)*0.0+...%dFSLHTM.ILM(j+1,i+1)=0.0 in Di-An system
                        RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1)+...
                        RS.PL(j+1,i+1)*HS.PL(j+1,i+1)*0.0;%dFSLHTM.PL(j+1,i+1)=0.0 when PL is not yet liquid phase, i.e., MCL.MgO>CE.MgO
                    
                    dFSLHN.CPX(j+1,i+1)=STLL/STLR;%[1]
                    dFSLHN.OL(j+1,i+1)=0.0;
                    dFSLHN.OPX(j+1,i+1)=0.0;
                    dFSLHN.PL(j+1,i+1)=0.0;
                    dFSLHN.ILM(j+1,i+1)=0.0;
                    
                    dFSLHN.OL(j+1,i+1)=0.5*(dFSLHN.OL(j+1,i+1)+dFSLHTM.OL(j+1,i+1));
                    dFSLHN.OPX(j+1,i+1)=0.5*(dFSLHN.OPX(j+1,i+1)+dFSLHTM.OPX(j+1,i+1));
                    dFSLHN.CPX(j+1,i+1)=0.5*(dFSLHN.CPX(j+1,i+1)+dFSLHTM.CPX(j+1,i+1));
                    dFSLHN.PL(j+1,i+1)=0.5*(dFSLHN.PL(j+1,i+1)+dFSLHTM.PL(j+1,i+1));
                    dFSLHN.ILM(j+1,i+1)=0.5*(dFSLHN.ILM(j+1,i+1)+dFSLHTM.ILM(j+1,i+1));
                    
                    errdFS(j+1,i+1)=max([abs((dFSLHN.OL(j+1,i+1)-dFSLHTM.OL(j+1,i+1))/dFSLHN.OL(j+1,i+1)),...
                        abs((dFSLHN.OPX(j+1,i+1)-dFSLHTM.OPX(j+1,i+1))/dFSLHN.OPX(j+1,i+1)),...
                        abs((dFSLHN.CPX(j+1,i+1)-dFSLHTM.CPX(j+1,i+1))/dFSLHN.CPX(j+1,i+1)),...
                        abs((dFSLHN.PL(j+1,i+1)-dFSLHTM.PL(j+1,i+1))/dFSLHN.PL(j+1,i+1)),...
                        abs((dFSLHN.ILM(j+1,i+1)-dFSLHTM.ILM(j+1,i+1))/dFSLHN.ILM(j+1,i+1))]);
                    
                    dFSLHTM.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1);
                    dFSLHTM.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1);
                    dFSLHTM.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1);
                    dFSLHTM.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1);
                    dFSLHTM.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1);
                    
                    %-------------------------- STEP THREE --------------------------------
                    %                     figure(2)
                    %                     plot(m,TN(j+1,i+1),'g*');
                    %                     hold on
                    %
                    %                     plot(m,TLN(j+1,i+1),'r*');
                    %                     hold on
                    %                     figure(3)
                    %                     plot(m,CLN(j+1,i+1),'k*');
                    %                     hold on
                    %                     figure(4)
                    %                     plot(m,dFSN(j+1,i+1),'ko')
                    %                     hold on
                    %                     m=m+1;
                    
                    %--------------------------- STEP FOUR --------------------------------
                    %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
                    %values, dFSLHN and CLN, and take TN=0.5*(TN+TTM).
                    
                    %NOTE: new T, CL have effects on RL, RS, CpS, CpL. For the first guess of T(n+1) and CL(n+1), which are differen from T(n) and CL(n),
                    %should drive updated RL, RS, CpS, CpL immediately. However, according to Temperature Eq., all but dFSLH are old values,
                    %so old RL, RS, CpS, CpS are used, and we don't update them here!!!
                    
                    %latent heat released in action [J/m^3]
                    LH.OL(j+1,i+1)=dFSLHN.OL(j+1,i+1)*RS.OL(j+1,i+1)*HS.OL(j+1,i+1);
                    LH.OPX(j+1,i+1)=dFSLHN.OPX(j+1,i+1)*RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1);
                    LH.CPX(j+1,i+1)=dFSLHN.CPX(j+1,i+1)*RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1);
                    LH.PL(j+1,i+1)=dFSLHN.PL(j+1,i+1)*RS.PL(j+1,i+1)*HS.PL(j+1,i+1);
                    LH.ILM(j+1,i+1)=dFSLHN.ILM(j+1,i+1)*RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1);
                    LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat release in action [J/m^3]
                    
                    %Old RS, CpS, New dFSLHN
                    RSCPSdFSTM(j+1,i+1)=RS.OL(j+1,i+1)*CPS.OL(j+1,i+1)*dFSLHN.OL(j+1,i+1)+...
                        RS.OPX(j+1,i+1)*CPS.OPX(j+1,i+1)*dFSLHN.OPX(j+1,i+1)+...
                        RS.CPX(j+1,i+1)*CPS.CPX(j+1,i+1)*dFSLHN.CPX(j+1,i+1)+...
                        RS.PL(j+1,i+1)*CPS.PL(j+1,i+1)*dFSLHN.PL(j+1,i+1)+...
                        RS.ILM(j+1,i+1)*CPS.ILM(j+1,i+1)*dFSLHN.ILM(j+1,i+1);
                    
                    dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                    TFN(j+1,i+1)=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1);%[J/m^3/K]
                    
                    TN(j+1,i+1)=(IHS(j+1,i+1)+TQVT(j,i)+TQDT(j,i)+LHT(j+1,i+1))/TFN(j+1,i+1);%[J/m^3 divided by J/m^3/K == K]
                    TN(j+1,i+1)=0.5*(TN(j+1,i+1)+TTM(j+1,i+1));
                    
                    errTTM(j+1,i+1)=abs((TN(j+1,i+1)-TTM(j+1,i+1))/TN(j+1,i+1));
                    
                    TTM(j+1,i+1)=TN(j+1,i+1);
                    %--------------------------- STEP FOUR --------------------------------
                    
                    %--------------------------- STEP FIVE --------------------------------
                    %With the new initial value TN (as well as dFSLHN and CLN), calculate the new approximation
                    %of CLN(n+1) using Eq. (9), and take CLN=0.5*(CLN+CLTM).
                    
                    %FL of of Next step
                    FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);%dFSSMT(j,i) probably ->0.0 since it's going to eutectic point
                    
                    %NOTE: new T has effect on RL,RS
                    %New RL (=RLN) updated by new T, MCL
                    [~,~,RLN(j+1,i+1)]=VisCpRL(TN(j+1,i+1),PATM(j+1,i+1),MCL0);
                    
                    %New CS updated by new T; originally New KP and then New CS in Xudaming
                    [~]=SYSEOS('Di');%MCSN
                    
                    %New RS (=RSN) updated by new T
                    RSNS=SolidDensity(TN(j+1,i+1),PATM(j+1,i+1));%solid phase density [kg/m^3]
                    RSN.OL(j+1,i+1)=RSNS.OL;
                    RSN.OPX(j+1,i+1)=RSNS.OPX;
                    RSN.CPX(j+1,i+1)=RSNS.CPX;
                    RSN.PL(j+1,i+1)=RSNS.PL;
                    RSN.ILM(j+1,i+1)=RSNS.ILM;
                    
                    %CLFNS--Factor ahead CL(n+1) for solid part: RSN, MCSN are updated of old RS, MCS;
                    %It should be:
                    %RSN.Fo(j+1,i+1)*NCSN.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+RS.Fo(j+1,i+1)*NCS.SiO2.Fo(j+1,i+1)*0.5*dFSLHTM.Fo(j+1,i+1)+...
                    %Try TWM, LPUM by pMELTS
                    
                    CLFNS.SiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.SiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.SiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.SiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.SiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.SiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.TiO2(j,i)=RSN.OL(j+1,i+1)*MCSN.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.TiO2.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.TiO2.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.TiO2.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.TiO2.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.TiO2.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Al2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Al2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Al2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Al2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Al2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Al2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.FeO(j,i)=RSN.OL(j+1,i+1)*MCSN.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.FeO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.FeO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.FeO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.FeO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.FeO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Fe2O3(j,i)=RSN.OL(j+1,i+1)*MCSN.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Fe2O3.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Fe2O3.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Fe2O3.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Fe2O3.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Fe2O3.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.MnO(j,i)=RSN.OL(j+1,i+1)*MCSN.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MnO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MnO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MnO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MnO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MnO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.MgO(j,i)=RSN.OL(j+1,i+1)*MCSN.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.MgO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.MgO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.MgO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.MgO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.MgO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.CaO(j,i)=RSN.OL(j+1,i+1)*MCSN.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.CaO.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.CaO.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.CaO.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.CaO.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.CaO.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Na2O(j,i)=RSN.OL(j+1,i+1)*MCSN.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.Na2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.Na2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.Na2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.Na2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.Na2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.K2O(j,i)=RSN.OL(j+1,i+1)*MCSN.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.K2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.K2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.K2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.K2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.K2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.P2O5(j,i)=RSN.OL(j+1,i+1)*MCSN.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.P2O5.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.P2O5.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.P2O5.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.P2O5.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.P2O5.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.H2O(j,i)=RSN.OL(j+1,i+1)*MCSN.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*MCS.H2O.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*MCSN.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*MCS.H2O.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*MCSN.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*MCS.H2O.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*MCSN.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*MCS.H2O.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*MCSN.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*MCS.H2O.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    %NOTE: For Trace elements, we use TCSN==TCS as approximation.
                    CLFNS.Sm(j,i)=RSN.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Sm.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Sm.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Sm.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Sm.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Sm.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    CLFNS.Nd(j,i)=RSN.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+RS.OL(j+1,i+1)*TCS.Nd.OL(j+1,i+1)*0.5*dFSLHTM.OL(j+1,i+1)+...
                        RSN.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+RS.OPX(j+1,i+1)*TCS.Nd.OPX(j+1,i+1)*0.5*dFSLHTM.OPX(j+1,i+1)+...
                        RSN.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+RS.CPX(j+1,i+1)*TCS.Nd.CPX(j+1,i+1)*0.5*dFSLHTM.CPX(j+1,i+1)+...
                        RSN.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+RS.PL(j+1,i+1)*TCS.Nd.PL(j+1,i+1)*0.5*dFSLHTM.PL(j+1,i+1)+...
                        RSN.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1)+RS.ILM(j+1,i+1)*TCS.Nd.ILM(j+1,i+1)*0.5*dFSLHTM.ILM(j+1,i+1);%[kg/m^3]
                    
                    %CLFN--Factor ahead CL(n+1) of liquid and solid: FL, RL should use NEW values FLN, RLN
                    CLFN.SiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.TiO2(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Al2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.FeO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Fe2O3(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.MnO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.MgO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.CaO(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Na2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.K2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.P2O5(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.H2O(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Sm(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    CLFN.Nd(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1);%[kg/m^3]
                    
                    %New major element CL of Next step in liquid [kg/m^3 divided by kg/m^3 ==1]
                    CLN.SiO2(j+1,i+1)=(MISS.SiO2(j+1,i+1)+CLQVT.SiO2(j,i)+CLQDT.SiO2(j,i)-CLFNS.SiO2(j,i))/CLFN.SiO2(j,i);
                    CLN.TiO2(j+1,i+1)=(MISS.TiO2(j+1,i+1)+CLQVT.TiO2(j,i)+CLQDT.TiO2(j,i)-CLFNS.TiO2(j,i))/CLFN.TiO2(j,i);
                    CLN.Al2O3(j+1,i+1)=(MISS.Al2O3(j+1,i+1)+CLQVT.Al2O3(j,i)+CLQDT.Al2O3(j,i)-CLFNS.Al2O3(j,i))/CLFN.Al2O3(j,i);
                    CLN.FeO(j+1,i+1)=(MISS.FeO(j+1,i+1)+CLQVT.FeO(j,i)+CLQDT.FeO(j,i)-CLFNS.FeO(j,i))/CLFN.FeO(j,i);
                    CLN.Fe2O3(j+1,i+1)=(MISS.Fe2O3(j+1,i+1)+CLQVT.Fe2O3(j,i)+CLQDT.Fe2O3(j,i)-CLFNS.Fe2O3(j,i))/CLFN.Fe2O3(j,i);
                    CLN.MnO(j+1,i+1)=(MISS.MnO(j+1,i+1)+CLQVT.MnO(j,i)+CLQDT.MnO(j,i)-CLFNS.MnO(j,i))/CLFN.MnO(j,i);
                    CLN.MgO(j+1,i+1)=(MISS.MgO(j+1,i+1)+CLQVT.MgO(j,i)+CLQDT.MgO(j,i)-CLFNS.MgO(j,i))/CLFN.MgO(j,i);
                    CLN.CaO(j+1,i+1)=(MISS.CaO(j+1,i+1)+CLQVT.CaO(j,i)+CLQDT.CaO(j,i)-CLFNS.CaO(j,i))/CLFN.CaO(j,i);
                    CLN.Na2O(j+1,i+1)=(MISS.Na2O(j+1,i+1)+CLQVT.Na2O(j,i)+CLQDT.Na2O(j,i)-CLFNS.Na2O(j,i))/CLFN.Na2O(j,i);
                    CLN.K2O(j+1,i+1)=(MISS.K2O(j+1,i+1)+CLQVT.K2O(j,i)+CLQDT.K2O(j,i)-CLFNS.K2O(j,i))/CLFN.K2O(j,i);
                    CLN.P2O5(j+1,i+1)=(MISS.P2O5(j+1,i+1)+CLQVT.P2O5(j,i)+CLQDT.P2O5(j,i)-CLFNS.P2O5(j,i))/CLFN.P2O5(j,i);
                    CLN.H2O(j+1,i+1)=(MISS.H2O(j+1,i+1)+CLQVT.H2O(j,i)+CLQDT.H2O(j,i)-CLFNS.H2O(j,i))/CLFN.H2O(j,i);
                    CLN.Sm(j+1,i+1)=(TISS.Sm(j+1,i+1)+CLQVT.Sm(j,i)+CLQDT.Sm(j,i)-CLFNS.Sm(j,i))/CLFN.Sm(j,i);
                    CLN.Nd(j+1,i+1)=(TISS.Nd(j+1,i+1)+CLQVT.Nd(j,i)+CLQDT.Nd(j,i)-CLFNS.Nd(j,i))/CLFN.Nd(j,i);
                    
                    %Take average CLN and CLTM
                    CLN.SiO2(j+1,i+1)=0.5*(CLN.SiO2(j+1,i+1)+CLTM.SiO2(j+1,i+1));
                    CLN.TiO2(j+1,i+1)=0.5*(CLN.TiO2(j+1,i+1)+CLTM.TiO2(j+1,i+1));
                    CLN.Al2O3(j+1,i+1)=0.5*(CLN.Al2O3(j+1,i+1)+CLTM.Al2O3(j+1,i+1));
                    CLN.FeO(j+1,i+1)=0.5*(CLN.FeO(j+1,i+1)+CLTM.FeO(j+1,i+1));
                    CLN.Fe2O3(j+1,i+1)=0.5*(CLN.Fe2O3(j+1,i+1)+CLTM.Fe2O3(j+1,i+1));
                    CLN.MnO(j+1,i+1)=0.5*(CLN.MnO(j+1,i+1)+CLTM.MnO(j+1,i+1));
                    CLN.MgO(j+1,i+1)=0.5*(CLN.MgO(j+1,i+1)+CLTM.MgO(j+1,i+1));
                    CLN.CaO(j+1,i+1)=0.5*(CLN.CaO(j+1,i+1)+CLTM.CaO(j+1,i+1));
                    CLN.Na2O(j+1,i+1)=0.5*(CLN.Na2O(j+1,i+1)+CLTM.Na2O(j+1,i+1));
                    CLN.K2O(j+1,i+1)=0.5*(CLN.K2O(j+1,i+1)+CLTM.K2O(j+1,i+1));
                    CLN.P2O5(j+1,i+1)=0.5*(CLN.P2O5(j+1,i+1)+CLTM.P2O5(j+1,i+1));
                    CLN.H2O(j+1,i+1)=0.5*(CLN.H2O(j+1,i+1)+CLTM.H2O(j+1,i+1));
                    CLN.Sm(j+1,i+1)=0.5*(CLN.Sm(j+1,i+1)+CLTM.Sm(j+1,i+1));
                    CLN.Nd(j+1,i+1)=0.5*(CLN.Nd(j+1,i+1)+CLTM.Nd(j+1,i+1));
                    
                    %Relative error of major elements in liquid [1]
                    errCLTM(j+1,i+1)=max([abs((CLN.SiO2(j+1,i+1)-CLTM.SiO2(j+1,i+1))/CLN.SiO2(j+1,i+1)),...
                        abs((CLN.TiO2(j+1,i+1)-CLTM.TiO2(j+1,i+1))/CLN.TiO2(j+1,i+1)),...
                        abs((CLN.Al2O3(j+1,i+1)-CLTM.Al2O3(j+1,i+1))/CLN.Al2O3(j+1,i+1)),...
                        abs((CLN.FeO(j+1,i+1)-CLTM.FeO(j+1,i+1))/CLN.FeO(j+1,i+1)),...
                        abs((CLN.Fe2O3(j+1,i+1)-CLTM.Fe2O3(j+1,i+1))/CLN.Fe2O3(j+1,i+1)),...
                        abs((CLN.MnO(j+1,i+1)-CLTM.MnO(j+1,i+1))/CLN.MnO(j+1,i+1)),...
                        abs((CLN.MgO(j+1,i+1)-CLTM.MgO(j+1,i+1))/CLN.MgO(j+1,i+1)),...
                        abs((CLN.CaO(j+1,i+1)-CLTM.CaO(j+1,i+1))/CLN.CaO(j+1,i+1)),...
                        abs((CLN.Na2O(j+1,i+1)-CLTM.Na2O(j+1,i+1))/CLN.Na2O(j+1,i+1)),...
                        abs((CLN.K2O(j+1,i+1)-CLTM.K2O(j+1,i+1))/CLN.K2O(j+1,i+1)),...
                        abs((CLN.P2O5(j+1,i+1)-CLTM.P2O5(j+1,i+1))/CLN.P2O5(j+1,i+1)),...
                        abs((CLN.H2O(j+1,i+1)-CLTM.H2O(j+1,i+1))/CLN.H2O(j+1,i+1))]);
                    
                    %Exchange old and new CL
                    CLTM.SiO2(j+1,i+1)=CLN.SiO2(j+1,i+1);
                    CLTM.TiO2(j+1,i+1)=CLN.TiO2(j+1,i+1);
                    CLTM.Al2O3(j+1,i+1)=CLN.Al2O3(j+1,i+1);
                    CLTM.FeO(j+1,i+1)=CLN.FeO(j+1,i+1);
                    CLTM.Fe2O3(j+1,i+1)=CLN.Fe2O3(j+1,i+1);
                    CLTM.MnO(j+1,i+1)=CLN.MnO(j+1,i+1);
                    CLTM.MgO(j+1,i+1)=CLN.MgO(j+1,i+1);
                    CLTM.CaO(j+1,i+1)=CLN.CaO(j+1,i+1);
                    CLTM.Na2O(j+1,i+1)=CLN.Na2O(j+1,i+1);
                    CLTM.K2O(j+1,i+1)=CLN.K2O(j+1,i+1);
                    CLTM.P2O5(j+1,i+1)=CLN.P2O5(j+1,i+1);
                    CLTM.H2O(j+1,i+1)=CLN.H2O(j+1,i+1);
                    CLTM.Sm(j+1,i+1)=CLN.Sm(j+1,i+1);
                    CLTM.Nd(j+1,i+1)=CLN.Nd(j+1,i+1);
                    
                    %--------------------------- STEP FIVE --------------------------------
                    
                    %Prepare Temporary MCL for VisCpRL.m
                    MCL0.SiO2=CLN.SiO2(j+1,i+1);
                    MCL0.TiO2=CLN.TiO2(j+1,i+1);
                    MCL0.Al2O3=CLN.Al2O3(j+1,i+1);
                    MCL0.FeO=CLN.FeO(j+1,i+1);
                    MCL0.Fe2O3=CLN.Fe2O3(j+1,i+1);
                    MCL0.MnO=CLN.MnO(j+1,i+1);
                    MCL0.MgO=CLN.MgO(j+1,i+1);
                    MCL0.CaO=CLN.CaO(j+1,i+1);
                    MCL0.Na2O=CLN.Na2O(j+1,i+1);
                    MCL0.K2O=CLN.K2O(j+1,i+1);
                    MCL0.P2O5=CLN.P2O5(j+1,i+1);
                    MCL0.H2O=CLN.H2O(j+1,i+1);
                    
                    %--------------------------- STEP SIX ---------------------------------
                    err(j+1,i+1)=logical((err0dFSLH<=errdFS(j+1,i+1))||(err0TTM<=errTTM(j+1,i+1))||(err0CLTM<=errCLTM(j+1,i+1)));%logical variable for while loop
                    
                    %If the conditions are not true, take TLN(n+1)=TLN(CLN(n+1)) and repeat steps 3-6 until
                    %all the above conditions are satisfied
                    
                    %--------------------------- STEP SIX ---------------------------------
                    
                    %--------------------------- STEP SEVEN -------------------------------
                    %If jugments are not true, then repeat step 3-6 with updated TN, CLN, dFSLHN
                    
                    %Critical major element taht represents system evolution
                    CLCM(j+1,i+1)=CLN.MgO(j+1,i+1);
                    
                    %Liquidus is a complex function of liquid composition, TLN=func(MgO, CaO, Al2O3, SiO2, H2O, ...)
                    TLN(j+1,i+1)=-1.3351*(CLCM(j+1,i+1)*100.0)^2+55.878*CLCM(j+1,i+1)*100.0+813.9991+273.15;%In Di-An system, new liquidus is a simple function of MgO [K]
                    
                    %-------------------------- STEP SEVEN --------------------------------
                    iters(j+1,i+1)=iters(j+1,i+1)+1;
                    
                end
                
                if(TLN(j+1,i+1)>TE)%back to case 2
                    FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);%dFSSMT(j,i) probably ->0.0 since it's going to eutectic point
                    %dFSLHN.ANY(j+1,i+1) are updated in while loop
                    %CLN.ANY(j+1,i+1) are updated in while loop
                    
                    %latent heat for energy check, updated in while loop [J/m^3]
                    LAT(j,i)=LHT(j+1,i+1);
                    
                    %sensible heat for energy check [J/m^3]
                    SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+(RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHT(j+1,i+1))*TN(j+1,i+1);
                    
                    %T for energy check [K]
                    TEC(j+1,i+1)=TN(j+1,i+1);
                    
                    %THERMAL EQUILIBRIUM --> T==Tliq
                    TN(j+1,i+1)=TLN(j+1,i+1);
                    
                    caseID(j,i)=2;%mark as situation 2
                    
                else%otherwise, case 3 is then valid
                    
                    %1. eutectic part solidification
                    QAE=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TTMBE(j+1,i+1)-TTMAE(j+1,i+1))*(TE-TLN(j+1,i+1))/(TTMBE(j+1,i+1)-TLN(j+1,i+1));%heat after eutectic
                    
                    %dFSLHNPLAE=0.599*RSE.PL*HS.CPX(j+1,i+1)/0.401+RSE.PL*HS.PL(j+1,i+1);
                    %0.599*RSE.PL*HS.CPX(j+1,i+1)/0.401+RSE.PL*HS.PL(j+1,i+1)+...%latent heat part
                    %0.599*RS.CPX(j+1,i+1)*CPS.CPX(j+1,i+1)*RSE.PL*TE/(0.401*RSE.CPX)+RS.PL(j+1,i+1)*CPS.PL(j+1,i+1)*TE-...%solid property change part
                    %(1.0+0.599*RSE.PL*RLTM(j+1,i+1)*CPLTM(j+1,i+1)/(0.401*RSE.CPX))*TE+...%solid property change part
                    %0.599*CPS.CPX(j+1,i+1)*RS.CPX(j+1,i+1)*RSE.PL*(TTMBE(j+1,i+1)-TE)/(0.401*RSE.CPX)+CPS.PL(j+1,i+1)*RS.PL(j+1,i+1)*(TTMBE(j+1,i+1)-TE);%sensible heat due to temperature change (TTME --> TE)
                    %NOTE: No sensible heat due to temperature change
                    
                    %new dFS of PL after eutectic [1]
                    dFSLHNPLAE=QAE/(0.599*RSE.PL*HS.CPX(j+1,i+1)/0.401+RSE.PL*HS.PL(j+1,i+1));%only latent heat
                    dFSLHNOLAE=0.0;
                    dFSLHNOPXAE=0.0;
                    dFSLHNILMAE=0.0;
                    dFSLHNCPXAE=0.599*dFSLHNPLAE*RSE.PL/(0.401*RSE.CPX);%new dFS of CPX after eutectic [1]
                    
                    %2. CPX cooling towards eutectic
                    QBE=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TTMBE(j+1,i+1)-TTMAE(j+1,i+1))*(TTMBE(j+1,i+1)-TE)/(TTMBE(j+1,i+1)-TLN(j+1,i+1));%heat before eutectic
                    %                     dFSLHNCPXBE=QBE/(RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1)+...%latent heat part
                    %                         (RS.CPX(j+1,i+1)*CPS.CPX(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1))*TE+...%sensible heat due to solid property change
                    %                         (FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TTMBE(j+1,i+1)-TE));%sensible heat due to temperature change
                    dFSLHNCPXBE=(QBE-(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TTMBE(j+1,i+1)-TE))/(RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1)-RS.CPX(j+1,i+1)*CPS.CPX(j+1,i+1)*TE+RLTM(j+1,i+1)*CPLTM(j+1,i+1)*TE);
                    dFSLHNOLBE=0.0;
                    dFSLHNPLBE=0.0;
                    dFSLHNOPXBE=0.0;
                    dFSLHNILMBE=0.0;
                    
                    dFSLHN.CPX(j+1,i+1)=dFSLHNCPXBE+dFSLHNCPXAE;
                    dFSLHN.PL(j+1,i+1)=dFSLHNPLBE+dFSLHNPLAE;
                    dFSLHN.OL(j+1,i+1)=dFSLHNOLBE+dFSLHNOLAE;
                    dFSLHN.OPX(j+1,i+1)=dFSLHNOPXBE+dFSLHNOPXAE;
                    dFSLHN.ILM(j+1,i+1)=dFSLHNILMBE+dFSLHNILMAE;
                    dFSLHT(j+1,i+1)=dFSLHN.OL(j+1,i+1)+dFSLHN.OPX(j+1,i+1)+dFSLHN.CPX(j+1,i+1)+dFSLHN.PL(j+1,i+1)+dFSLHN.ILM(j+1,i+1);
                    
                    dFSLHNBE=dFSLHNCPXBE+dFSLHNOPXBE+dFSLHNOLBE+dFSLHNPLBE+dFSLHNILMBE;%total dFSLH prior to eutectic
                    dFSLHNAE=dFSLHNCPXAE+dFSLHNOPXAE+dFSLHNOLAE+dFSLHNPLAE+dFSLHNILMAE;%total dFSLH after to eutectic
                    
                    %latent heat released in action [J/m^3]
                    LH.OL(j+1,i+1)=dFSLHNOLBE*RS.OL(j+1,i+1)*HS.OL(j+1,i+1)+...%TTMBE --> TE latent heat
                        dFSLHNOLAE*RSE.OL*HS.OL(j+1,i+1);%TE latent heat
                    LH.OPX(j+1,i+1)=dFSLHNOPXBE*RS.OPX(j+1,i+1)*HS.OPX(j+1,i+1)+...%TTMBE --> TE latent heat
                        dFSLHNOPXAE*RSE.OPX*HS.OPX(j+1,i+1);%TE latent heat
                    LH.CPX(j+1,i+1)=dFSLHNCPXBE*RS.CPX(j+1,i+1)*HS.CPX(j+1,i+1)+...%TTMBE --> TE latent heat
                        dFSLHNCPXAE*RSE.CPX*HS.CPX(j+1,i+1);%TE latent heat
                    LH.PL(j+1,i+1)=dFSLHNPLBE*RS.PL(j+1,i+1)*HS.PL(j+1,i+1)+...%TTMBE --> TE latent heat
                        dFSLHNPLAE*RSE.PL*HS.PL(j+1,i+1);%TE latent heat
                    LH.ILM(j+1,i+1)=dFSLHNILMBE*RS.ILM(j+1,i+1)*HS.ILM(j+1,i+1)+...%TTMBE --> TE latent heat
                        dFSLHNILMAE*RSE.ILM*HS.ILM(j+1,i+1);%TE latent heat
                    
                    LHT(j+1,i+1)=LH.OL(j+1,i+1)+LH.OPX(j+1,i+1)+LH.CPX(j+1,i+1)+LH.PL(j+1,i+1)+LH.ILM(j+1,i+1);%total latent heat release in action [J/m^3]
                    
                    %dFSSEN(j+1,i+1)=dFSCAL;%used in energy check for sensible heat
                    modifymarker=modifymarker+1;
                    
                    FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSLHT(j+1,i+1);%dFSSMT(j,i) probably ->0.0 since it's going to eutectic point
                    
                    TN(j+1,i+1)=TE;
                    
                    CLN.SiO2(j+1,i+1)=CE.SiO2;
                    CLN.TiO2(j+1,i+1)=CE.TiO2;
                    CLN.Al2O3(j+1,i+1)=CE.Al2O3;
                    CLN.FeO(j+1,i+1)=CE.FeO;
                    CLN.Fe2O3(j+1,i+1)=CE.Fe2O3;
                    CLN.MnO(j+1,i+1)=CE.MnO;
                    CLN.MgO(j+1,i+1)=CE.MgO;
                    CLN.CaO(j+1,i+1)=CE.CaO;
                    CLN.Na2O(j+1,i+1)=CE.Na2O;
                    CLN.K2O(j+1,i+1)=CE.K2O;
                    CLN.P2O5(j+1,i+1)=CE.P2O5;
                    CLN.H2O(j+1,i+1)=CE.H2O;
                    %NOTE: We use last step TCL as eutectic TCL for simplicity.
                    CLN.Sm(j+1,i+1)=TCL.Sm(j+1,i+1);
                    CLN.Nd(j+1,i+1)=TCL.Nd(j+1,i+1);
                    
                    %latent heat for energy check [J/m^3]
                    LAT(j,i)=LHT(j+1,i+1);
                    
                    %[1] sensible heat before eutectic for energy check [J/m^3]
                    %Old RS, CpS, New dFSLHN
                    RSCPSdFSTM(j+1,i+1)=RS.OL(j+1,i+1)*CPS.OL(j+1,i+1)*dFSLHNOLBE+...
                        RS.OPX(j+1,i+1)*CPS.OPX(j+1,i+1)*dFSLHNOPXBE+...
                        RS.CPX(j+1,i+1)*CPS.CPX(j+1,i+1)*dFSLHNCPXBE+...
                        RS.PL(j+1,i+1)*CPS.PL(j+1,i+1)*dFSLHNPLBE+...
                        RS.ILM(j+1,i+1)*CPS.ILM(j+1,i+1)*dFSLHNILMBE;
                    %NOTE: while loop in case 3 doesn't give physically right answer, but dFSLHN___BE and dFSLHN___AE are right
                    
                    SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TN(j+1,i+1)-TTMBE(j+1,i+1))+...%sensible heat due to temperature change during TTME --> TE; TN=TE
                        (RSCPSdFSTM(j+1,i+1)-RLTM(j+1,i+1)*CPLTM(j+1,i+1)*dFSLHNBE)*TN(j+1,i+1);%sensible heat due to solid property change during TTME --> TE
                    
                    %                     %[2] sensible heat after eutectic for energy check [J/m^3]
                    %                     RSCPSdFSTM(j+1,i+1)=RSE.OL*CPSE.OL*dFSLHNOLAE+...
                    %                         RSE.OPX*CPSE.OPX*dFSLHNOPXAE+...
                    %                         RSE.CPX*CPSE.CPX*dFSLHNCPXAE+...
                    %                         RSE.PL*CPSE.PL*dFSLHNPLAE+...
                    %                         RSE.ILM*CPSE.ILM*dFSLHNILMAE;
                    %                     SEN(j,i)=SEN(j,i)+0.0+(RSCPSdFSTM(j+1,i+1)-RLE*CPLE*dFSLHNAE)*TN(j+1,i+1);
                    %NOTE: TN==TE, '0.0' is sensible heat due to temperature change
                    
                    %T for energy check [K]
                    TEC(j+1,i+1)=TN(j+1,i+1);
                    
                    if(FSELOG(j+1,i+1))%FSE only modify once
                        FSE(j+1,i+1)=FSTTM(j+1,i+1)+dFSLHNBE;%real solid fraction when eutectic
                        FSELOG(j+1,i+1)=-1.0;
                    end
                    
                    fprintf('=============== FSE(%2d,%2d)= %7.5f ==============\n',j+1,i+1,FSE(j+1,i+1));
                end
                
                
                %                 else
                %                     fprintf('Complete solid BUT should be eutectic\n');
                %                 end
                
        end
        
    end
end
%################################################################ T-FS-CL ITERATION #########################################################################


%% ========================== ENERGY CHECK ================================

%Absolute difference between energy removed and energy loss, RES=RESidual, E=Energy
RESE=zeros(NIY+2,NIX+2);%[J/m^3]
RESER=zeros(NIY,NIX);
%SENMAX=max(max(abs(SEN(1:NIY,1:NIX))));
energy_check=1.0;%energy check marker

for i=1:NIX
    for j=1:NIY
        RESE(j+1,i+1)=SEN(j,i)-TQVT(j,i)-LAT(j,i)-TQDT(j,i);%All in [J/m^3]
        %NOTE: RESE does not include geometry configuration, i.e., dx, dy
        %NOTE: better use absolute value
        RESER(j,i)=abs(RESE(j+1,i+1)/SEN(j,i));
    end
end

RESE(1,2:NIX+1)=RESE(2,2:NIX+1);
RESE(NIY+2,2:NIX+1)=RESE(NIY+1,2:NIX+1);
RESE(1:NIY+2,1)=RESE(1:NIY+2,2);
RESE(1:NIY+2,NIX+2)=RESE(1:NIY+2,NIX+1);
MaxErr=zeros(NIY,NIX);

%Scheme 1: using significant figures to set error limit
for i=1:NIX
    for j=1:NIY
        if((abs(RESE(j+1,i+1))>=1.0e-5)&&(RESER(j,i)>=1.0e-5))
            fprintf('Energy Check: Cell (%2d,%2d) not balanced; AE: %E, RE: %E\n',j,i,RESER(j,i),RESE(j+1,i+1));
            energy_check=-1.0;
        end
        MaxErr(j,i)=min(abs(RESER(j,i)),abs(RESE(j+1,i+1)));
    end
end

if(energy_check>0.0)
    fprintf('Energy Balanced -- Max error %E\n',max(max(MaxErr)));
end

%Scheme 2: using relative error to set error limit
% for i=1:NIX
%     for j=1:NIY
%
%         if(abs(TQDT(j,i)/TQDTMAX)<=1.0e-12)%machine error of TQDT
%             TQDT(j,i)=0.0;
%         end
%
%         if(isinf(RESER(j,i)))%RESE~=0 but TQDT=0
%             if(abs(RESE(j+1,i+1))>=1.0e-6)%machine error but inevitable
%                 fprintf('Machine Error: %E at Cell (%2d,%2d)\n',RESE(j+1,i+1),j,i);
%             else
%                 RESER(j,i)=0.0;
%             end
%         end
%
%         if(isnan(RESER(j,i)))%RESE==0 and TQDT==0
%             RESER(j,i)=0.0;
%         end
%
%         if(RESER(j,i)>=1.0e-6)%NOTE: 1.0e-6 is a very coarse criteria. Most of them are below 1.0e-10
%             fprintf('Energy Error: %E at Cell (%2d,%2d)\n',RESER(j,i),j,i);
%             energy_check=-1.0;
%         end
%     end
% end
%
% if(energy_check>0.0)
%     fprintf('Energy Balanced -- Max error %E\n',max(max(RESER)));
% end
%---------------------------- ENERGY CHECK --------------------------------

%% ============== UPDATE MAJOR & TRACE ELEMENTS CONCENTRATION =============

%[1] Update Major/Trace Residual Internal Species Storage [kg/m^3]
for i=2:NIX+1
    for j=2:NIY+1
        %NOTE: CLFNS.ANY are updated from main iteration above and the last CLFNS.ANY meet the exit requirements and are thus considered as the results
        MISS.SiO2(j,i)=MISS.SiO2(j,i)+CLQVT.SiO2(j-1,i-1)+CLQDT.SiO2(j-1,i-1)-CLFNS.SiO2(j-1,i-1);%[kg/m^3]
        MISS.TiO2(j,i)=MISS.TiO2(j,i)+CLQVT.TiO2(j-1,i-1)+CLQDT.TiO2(j-1,i-1)-CLFNS.TiO2(j-1,i-1);%[kg/m^3]
        MISS.Al2O3(j,i)=MISS.Al2O3(j,i)+CLQVT.Al2O3(j-1,i-1)+CLQDT.Al2O3(j-1,i-1)-CLFNS.Al2O3(j-1,i-1);%[kg/m^3]
        MISS.FeO(j,i)=MISS.FeO(j,i)+CLQVT.FeO(j-1,i-1)+CLQDT.FeO(j-1,i-1)-CLFNS.FeO(j-1,i-1);%[kg/m^3]
        MISS.Fe2O3(j,i)=MISS.Fe2O3(j,i)+CLQVT.Fe2O3(j-1,i-1)+CLQDT.Fe2O3(j-1,i-1)-CLFNS.Fe2O3(j-1,i-1);%[kg/m^3]
        MISS.MnO(j,i)=MISS.MnO(j,i)+CLQVT.MnO(j-1,i-1)+CLQDT.MnO(j-1,i-1)-CLFNS.MnO(j-1,i-1);%[kg/m^3]
        MISS.MgO(j,i)=MISS.MgO(j,i)+CLQVT.MgO(j-1,i-1)+CLQDT.MgO(j-1,i-1)-CLFNS.MgO(j-1,i-1);%[kg/m^3]
        MISS.CaO(j,i)=MISS.CaO(j,i)+CLQVT.CaO(j-1,i-1)+CLQDT.CaO(j-1,i-1)-CLFNS.CaO(j-1,i-1);%[kg/m^3]
        MISS.Na2O(j,i)=MISS.Na2O(j,i)+CLQVT.Na2O(j-1,i-1)+CLQDT.Na2O(j-1,i-1)-CLFNS.Na2O(j-1,i-1);%[kg/m^3]
        MISS.K2O(j,i)=MISS.K2O(j,i)+CLQVT.K2O(j-1,i-1)+CLQDT.K2O(j-1,i-1)-CLFNS.K2O(j-1,i-1);%[kg/m^3]
        MISS.P2O5(j,i)=MISS.P2O5(j,i)+CLQVT.P2O5(j-1,i-1)+CLQDT.P2O5(j-1,i-1)-CLFNS.P2O5(j-1,i-1);%[kg/m^3]
        MISS.H2O(j,i)=MISS.H2O(j,i)+CLQVT.H2O(j-1,i-1)+CLQDT.H2O(j-1,i-1)-CLFNS.H2O(j-1,i-1);%[kg/m^3]
        
        TISS.Sm(j,i)=TISS.Sm(j,i)+CLQVT.Sm(j-1,i-1)+CLQDT.Sm(j-1,i-1)-CLFNS.Sm(j-1,i-1);%[kg/m^3]
        TISS.Nd(j,i)=TISS.Nd(j,i)+CLQVT.Nd(j-1,i-1)+CLQDT.Nd(j-1,i-1)-CLFNS.Nd(j-1,i-1);%[kg/m^3]
    end
end
MISS.SiO2=Border(MISS.SiO2);
MISS.TiO2=Border(MISS.TiO2);
MISS.Al2O3=Border(MISS.Al2O3);
MISS.FeO=Border(MISS.FeO);
MISS.Fe2O3=Border(MISS.Fe2O3);
MISS.MnO=Border(MISS.MnO);
MISS.MgO=Border(MISS.MgO);
MISS.CaO=Border(MISS.CaO);
MISS.Na2O=Border(MISS.Na2O);
MISS.K2O=Border(MISS.K2O);
MISS.P2O5=Border(MISS.P2O5);
MISS.H2O=Border(MISS.H2O);

TISS.Sm=Border(TISS.Sm);
TISS.Nd=Border(TISS.Nd);

%[2] Update dFSLH[NIY+2,NIX+2]
%dFSLHN.ANY(1,:) must be 0.0 to match real top cold boundary since RSCPSdFSTM, dFSLHT, TFN, LH.ANY and LHT depend on dFSLHN.ANY
%Recheck dFSLHN.ANY(1,:)=0.0
dFSLHN.OL(1,:)=0.0;
dFSLHN.OPX(1,:)=0.0;
dFSLHN.CPX(1,:)=0.0;
dFSLHN.PL(1,:)=0.0;
dFSLHN.ILM(1,:)=0.0;

dFSLH.OL=Border(dFSLHN.OL);
%NOTE: dFSLH.OL=Border(dFSLHN.OL) is equal to the following:
%     dFSLHN.OL(1,2:NIX+1)=0.0;%top cool boundary, stays the same
%     dFSLHN.OL(NIY+2,2:NIX+1)=dFSLHN.OL(NIY+1,2:NIX+1);%bottom insulated boundary
%     dFSLHN.OL(1:NIY+2,1)=dFSLHN.OL(1:NIY+2,2);%left insulated boundary
%     dFSLHN.OL(1:NIY+2,NIX+2)=dFSLHN.OL(1:NIY+2,NIX+1);%right insulated boundary
%     dFSLH.OL=dFSLHN.OL;

dFSLH.OPX=Border(dFSLHN.OPX);
dFSLH.CPX=Border(dFSLHN.CPX);
dFSLH.PL=Border(dFSLHN.PL);
dFSLH.ILM=Border(dFSLHN.ILM);

%[3] Update dFSSM
%dFSSM[NIY+2,NIX+2] and dFSSMT[NIY+2,NIX+2] are updated in SolidAssemble.m before MAGTFC.m

%[4] Update dFS[NIY+2,NIX+2]
dFS.OL=dFSLH.OL;
dFS.OPX=dFSLH.OPX;
dFS.CPX=dFSLH.CPX;
dFS.PL=dFSLH.PL;
dFS.ILM=dFSLH.ILM;

%[5] Update FS[NIY+2,NIX+2]
FS.OL=FS.OL+dFS.OL;
FS.OPX=FS.OPX+dFS.OPX;
FS.CPX=FS.CPX+dFS.CPX;
FS.PL=FS.PL+dFS.PL;
FS.ILM=FS.ILM+dFS.ILM;

%[6] Update MCS[NIY+2,NIX+2]
MCS.SiO2.OL=Border(MCSN.SiO2.OL);
MCS.SiO2.OPX=Border(MCSN.SiO2.OPX);
MCS.SiO2.CPX=Border(MCSN.SiO2.CPX);
MCS.SiO2.PL=Border(MCSN.SiO2.PL);
MCS.SiO2.ILM=Border(MCSN.SiO2.ILM);
%NOTE: MCL.SiO2=Border(CLN.SiO2) is equal to the following:
%     MCL.SiO2(1,2:NIX+1)=MCL.SiO2(1,2:NIX+1)_old;%top cool boundary, stays the same
%     MCL.SiO2(NIY+2,2:NIX+1)=MCL.SiO2(NIY+1,2:NIX+1);%bottom boundary
%     MCL.SiO2(1:NIY+2,1)=MCL.SiO2(1:NIY+2,2);%left boundary
%     MCL.SiO2(1:NIY+2,NIX+2)=MCL.SiO2(1:NIY+2,NIX+1);%right boundary

MCS.TiO2.OL=Border(MCSN.TiO2.OL);
MCS.TiO2.OPX=Border(MCSN.TiO2.OPX);
MCS.TiO2.CPX=Border(MCSN.TiO2.CPX);
MCS.TiO2.PL=Border(MCSN.TiO2.PL);
MCS.TiO2.ILM=Border(MCSN.TiO2.ILM);

MCS.Al2O3.OL=Border(MCSN.Al2O3.OL);
MCS.Al2O3.OPX=Border(MCSN.Al2O3.OPX);
MCS.Al2O3.CPX=Border(MCSN.Al2O3.CPX);
MCS.Al2O3.PL=Border(MCSN.Al2O3.PL);
MCS.Al2O3.ILM=Border(MCSN.Al2O3.ILM);

MCS.FeO.OL=Border(MCSN.FeO.OL);
MCS.FeO.OPX=Border(MCSN.FeO.OPX);
MCS.FeO.CPX=Border(MCSN.FeO.CPX);
MCS.FeO.PL=Border(MCSN.FeO.PL);
MCS.FeO.ILM=Border(MCSN.FeO.ILM);

MCS.Fe2O3.OL=Border(MCSN.Fe2O3.OL);
MCS.Fe2O3.OPX=Border(MCSN.Fe2O3.OPX);
MCS.Fe2O3.CPX=Border(MCSN.Fe2O3.CPX);
MCS.Fe2O3.PL=Border(MCSN.Fe2O3.PL);
MCS.Fe2O3.ILM=Border(MCSN.Fe2O3.ILM);

MCS.MnO.OL=Border(MCSN.MnO.OL);
MCS.MnO.OPX=Border(MCSN.MnO.OPX);
MCS.MnO.CPX=Border(MCSN.MnO.CPX);
MCS.MnO.PL=Border(MCSN.MnO.PL);
MCS.MnO.ILM=Border(MCSN.MnO.ILM);

MCS.MgO.OL=Border(MCSN.MgO.OL);
MCS.MgO.OPX=Border(MCSN.MgO.OPX);
MCS.MgO.CPX=Border(MCSN.MgO.CPX);
MCS.MgO.PL=Border(MCSN.MgO.PL);
MCS.MgO.ILM=Border(MCSN.MgO.ILM);

MCS.CaO.OL=Border(MCSN.CaO.OL);
MCS.CaO.OPX=Border(MCSN.CaO.OPX);
MCS.CaO.CPX=Border(MCSN.CaO.CPX);
MCS.CaO.PL=Border(MCSN.CaO.PL);
MCS.CaO.ILM=Border(MCSN.CaO.ILM);

MCS.Na2O.OL=Border(MCSN.Na2O.OL);
MCS.Na2O.OPX=Border(MCSN.Na2O.OPX);
MCS.Na2O.CPX=Border(MCSN.Na2O.CPX);
MCS.Na2O.PL=Border(MCSN.Na2O.PL);
MCS.Na2O.ILM=Border(MCSN.Na2O.ILM);

MCS.K2O.OL=Border(MCSN.K2O.OL);
MCS.K2O.OPX=Border(MCSN.K2O.OPX);
MCS.K2O.CPX=Border(MCSN.K2O.CPX);
MCS.K2O.PL=Border(MCSN.K2O.PL);
MCS.K2O.ILM=Border(MCSN.K2O.ILM);

MCS.P2O5.OL=Border(MCSN.P2O5.OL);
MCS.P2O5.OPX=Border(MCSN.P2O5.OPX);
MCS.P2O5.CPX=Border(MCSN.P2O5.CPX);
MCS.P2O5.PL=Border(MCSN.P2O5.PL);
MCS.P2O5.ILM=Border(MCSN.P2O5.ILM);

MCS.H2O.OL=Border(MCSN.H2O.OL);
MCS.H2O.OPX=Border(MCSN.H2O.OPX);
MCS.H2O.CPX=Border(MCSN.H2O.CPX);
MCS.H2O.PL=Border(MCSN.H2O.PL);
MCS.H2O.ILM=Border(MCSN.H2O.ILM);

%[7] Update RS[NIY+2,NIX+2]
RS.OL=Border(RSN.OL);
RS.OPX=Border(RSN.OPX);
RS.CPX=Border(RSN.CPX);
RS.PL=Border(RSN.PL);
RS.ILM=Border(RSN.ILM);

%[8] Update MCL[NIY+2,NIX+2]
MCL.SiO2=Border(CLN.SiO2);
MCL.TiO2=Border(CLN.TiO2);
MCL.Al2O3=Border(CLN.Al2O3);
MCL.FeO=Border(CLN.FeO);
MCL.Fe2O3=Border(CLN.Fe2O3);
MCL.MnO=Border(CLN.MnO);
MCL.MgO=Border(CLN.MgO);
MCL.CaO=Border(CLN.CaO);
MCL.Na2O=Border(CLN.Na2O);
MCL.K2O=Border(CLN.K2O);
MCL.P2O5=Border(CLN.P2O5);
MCL.H2O=Border(CLN.H2O);

MCL.SiO2(1,:)=MCL.SiO2(2,:);
MCL.TiO2(1,:)=MCL.TiO2(2,:);
MCL.Al2O3(1,:)=MCL.Al2O3(2,:);
MCL.FeO(1,:)=MCL.FeO(2,:);
MCL.Fe2O3(1,:)=MCL.Fe2O3(2,:);
MCL.MnO(1,:)=MCL.MnO(2,:);
MCL.MgO(1,:)=MCL.MgO(2,:);
MCL.CaO(1,:)=MCL.CaO(2,:);
MCL.Na2O(1,:)=MCL.Na2O(2,:);
MCL.K2O(1,:)=MCL.K2O(2,:);
MCL.P2O5(1,:)=MCL.P2O5(2,:);
MCL.H2O(1,:)=MCL.H2O(2,:);


%[9] Update TCL[NIY+2,NIX+2]
TCL.Sm=Border(CLN.Sm);
TCL.Nd=Border(CLN.Nd);
TCL.Sm(1,:)=TCL.Sm(2,:);
TCL.Nd(1,:)=TCL.Nd(2,:);

for i=1:NIX+2
    for j=1:NIY+2
        
        %[11] Update Total major elements in solid at each step (SM+LH) [NIY+2,NIX+2]
        %-------------------------- FORTRAN -------------------------------
        FRCSMj.SiO2(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.SiO2.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.SiO2.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.SiO2.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.SiO2.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.SiO2.ILM(j,i);
        FRCSMj.TiO2(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.TiO2.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.TiO2.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.TiO2.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.TiO2.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.TiO2.ILM(j,i);
        FRCSMj.Al2O3(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.Al2O3.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.Al2O3.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.Al2O3.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.Al2O3.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.Al2O3.ILM(j,i);
        FRCSMj.FeO(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.FeO.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.FeO.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.FeO.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.FeO.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.FeO.ILM(j,i);
        FRCSMj.Fe2O3(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.Fe2O3.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.Fe2O3.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.Fe2O3.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.Fe2O3.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.Fe2O3.ILM(j,i);
        FRCSMj.MnO(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.MnO.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.MnO.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.MnO.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.MnO.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.MnO.ILM(j,i);
        FRCSMj.MgO(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.MgO.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.MgO.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.MgO.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.MgO.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.MgO.ILM(j,i);
        FRCSMj.CaO(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.CaO.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.CaO.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.CaO.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.CaO.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.CaO.ILM(j,i);
        FRCSMj.Na2O(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.Na2O.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.Na2O.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.Na2O.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.Na2O.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.Na2O.ILM(j,i);
        FRCSMj.K2O(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.K2O.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.K2O.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.K2O.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.K2O.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.K2O.ILM(j,i);
        FRCSMj.P2O5(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.P2O5.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.P2O5.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.P2O5.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.P2O5.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.P2O5.ILM(j,i);
        FRCSMj.H2O(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*MCS.H2O.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*MCS.H2O.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*MCS.H2O.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*MCS.H2O.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*MCS.H2O.ILM(j,i);
        
        %[12] Update KSL[NIY+2,NIX+2]
        %KSL may change due to new T and new TCL and update here
        KP.Sm.CPX(j,i)=0.3;%maybe a function of T and TCL.Sm
        KP.Nd.CPX(j,i)=0.2;%maybe a function of T and TCL.Nd
        %...more KP ...
        KP.Sm.PL(j,i)=0.03;%maybe a function of T and TCL.Sm
        KP.Nd.PL(j,i)=0.04;%maybe a function of T and TCL.Nd
        
        
        %[13] Update TCS[NIY+2,NIX+2]
        TCS.Sm.OL(j,i)=TCL.Sm(j,i)*KP.Sm.OL(j,i);
        TCS.Nd.OL(j,i)=TCL.Nd(j,i)*KP.Nd.OL(j,i);
        %...more TCS ...
        TCS.Sm.PL(j,i)=TCL.Sm(j,i)*KP.Sm.OL(j,i);
        TCS.Nd.PL(j,i)=TCL.Nd(j,i)*KP.Nd.OL(j,i);
        
        
        %[14] Update Total trace elements in solid at each step (SM+LH) [NIY+2,NIX+2]
        FRCSTe.Sm(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*TCS.Sm.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*TCS.Sm.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*TCS.Sm.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*TCS.Sm.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*TCS.Sm.ILM(j,i);
        FRCSTe.Nd(j,i)=dFSLH.OL(j,i)*RS.OL(j,i)*TCS.Nd.OL(j,i)+dFSLH.OPX(j,i)*RS.OPX(j,i)*TCS.Nd.OPX(j,i)+dFSLH.CPX(j,i)*RS.CPX(j,i)*TCS.Nd.CPX(j,i)+dFSLH.PL(j,i)*RS.PL(j,i)*TCS.Nd.PL(j,i)+dFSLH.ILM(j,i)*RS.ILM(j,i)*TCS.Nd.ILM(j,i);
        
    end
end

%[15] Update all T, FL
%Returen new field variables
TN(1,2:NIX+1)=TW;%top real cool boundary
TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%bottom insulated boundary
TN(1:NIY+2,1)=TN(1:NIY+2,2);%left insulated boundary
TN(1:NIY+2,NIX+2)=TN(1:NIY+2,NIX+1);%right insulated boundary
T=TN;

mf0=80.0/216.5504/(80.0/216.5504+20.0/278.2073);%Di initial mole fraction
KSM=(mf0*KST.CPX+(1.0-mf0)*KST.PL)+0.0*KST.OL+0.0*KST.OPX+0.0*KST.ILM;
Q=KSM*(TN(2,i)-TW)/(0.5*dy(1));
QM=sum(Q)/NIX;%[W.m^-2]

FL=Border(FLN);
%Recheck real top boundary
FL(1,:)=0.0;
%NOTE: FL=Border(FLN) is equal to the following:
%     FLN(1,1:NIX+2)=0.0;%top cool boundary, stays the same
%     FLN(NIY+2,2:NIX+1)=FLN(NIY+1,2:NIX+1);%bottom boundary
%     FLN(1:NIY+2,1)=FLN(1:NIY+1,2);%left boundary
%     FLN(1:NIY+2,NIX+2)=FLN(1:NIY+2,NIX+1);%right boundary
%     FL=FLN;

fprintf('T-FS-CL Loop:%3d iterations!\n',max(max(iters)));

end