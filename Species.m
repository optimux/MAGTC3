function [CLQDT,CLQVT]=Species(CLTM,DLTM,FLTM,RLTM,RFVX,RFVY)
%This function is used to calculate species conservation of both major and trace elements
%CLQDT, CLQVT, CSQVT: output of old-time values in species conservation equation for all oxide, i.e., right hand side of species equation
%Oxide: name of oxide, input as a string
%Created on 2020-7-2

global NIX
global NIY
global dtb
global dx
global dy
global MCL
global TCL
global DL
global MCS
global TCS
%global ISS

%CLTM=zeros(NIY+2,NIX+2);%general CL
%DLTM=zeros(NIY+2,NIX+2);%general DL, and DL is updated during each time step by LiqDiff.m

%% =================== DIFFUSION SPECIES FLUX ====================

%x-axis FL*RL*CL
FRCLX=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        FRCLX(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
    end
end

%y-axis FL*RL*CL
FRCLY=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        FRCLY(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
    end
end

%x-axis diffusion
CLQDX=zeros(NIY+2,NIX+1);%CL=CL, Q=flux, D=diffusion, X=x-axis
%VERY IMPORTANT NOTE: CLQDX!=0.0 if CL(j,i+1)==CL(j,i) but RLTM(j,i+1)!=RLTM(j,i). This happens
%when no crystallization occurs but temperature keeps decreasing which led
%to increasing liquid density. However, diffusion flux of this is orders of
%magnitute smaller than other mechanisms, we therefore let this be here.
%BUT, if some species are highly mobile, we'd use rigorously concentration
%gradient!
for i=2:NIX
    %CL(1:NIY+2,1)=CL(1:NIY+2,2) --> CLQDX(1:NIY+2,1)=0.0 (left impermeable boundary)
    %CL(1:NIY+2,NIX+2)=CL(1:NIY+2,NIX+1) --> CLQDX(1:NIY+2,NIX+1)=0.0 (right impermeable boundary)
    for j=1:NIY+2
        CLQDX(j,i)=2.0*dtb*0.5*(DLTM(j,i)+DLTM(j,i+1))*(FRCLX(j,i+1)-FRCLX(j,i))/(dx(i-1)+dx(i));%[kg/m^2]
    end
end

%y-axis diffusion
CLQDY=zeros(NIY+1,NIX+2);%CL=CL, Q=flux, D=diffusion, Y=y-axis
%VERY IMPORTANT NOTE: CLQDY!=0.0 if CL(j,i)==CL(j+1,i) but RLTM(j,i)!=RLTM(j+1,i). This happens
%when no crystallization occurs but temperature keeps decreasing which led
%to increasing liquid density. However, diffusion flux of this is orders of
%magnitute smaller than other mechanisms, we therefore let this be here.
%BUT, if some species are highly mobile, we'd use rigorously concentration
%gradient!
for i=1:NIX+2
    for j=2:NIY
        %CL(1,1:NIX+2)=CL(2,1:NIX+2) --> CLQDY(1,1:NIX+2)=0.0 (top impermeable boundary)
        %CL(NIY+2,1:NIX+2)=CL(NIY+1,1:NIX+2) --> CLQDY(NIY+1,1:NIX+2)=0.0 (bottom impermeable boundary)
        CLQDY(j,i)=2.0*dtb*0.5*(DLTM(j,i)+DLTM(j+1,i))*(FRCLY(j+1,i)-FRCLY(j,i))/(dy(j-1)+dy(j));%[kg/m^2]
    end
end

%total species flux of each finite volume
CLQDT=zeros(NIY,NIX);%CL=CL, Q=flux, D=diffusion, T=total
for i=1:NIX
    for j=1:NIY
        CLQDT(j,i)=(CLQDX(j+1,i+1)-CLQDX(j+1,i))/dx(i)+(CLQDY(j+1,i+1)-CLQDY(j,i+1))/dy(j);%[kg/m^3]
    end
end

%% =============== LIQUID CONVECTION SPECIES FLUX =================
%---------------------- Part One: CLRFVX CLRFVY ---------------------------

%x-axis species liquid convection flux
CLRFVX=zeros(NIY+2,NIX+1);
for i=2:NIX
    %VX(1:NIY+2,1)=0.0 --> CLRFVX(1:NIY+2,1)=0.0 (left impermeable boundary)
    %VX(1:NIY+2,NIX+1)=0.0 --> CLRFVX(1:NIY+2,NIX+1)=0.0 (right impermeable boundary)
    for j=1:NIY+2
        %VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> CLRFVX(1,1:NIX+1)=-CLRFVX(2,1:NIX+1) --> No mass flux along a-axis boundary --> top NO SLIP boundary
        %VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1) --> CLRFVX(NIY+2,1:NIX+1)=-CLRFVX(NIY+1,1:NIX+1) --> No mass flux along x-axis boundary --> bottom NO SLIP boundary
        
        %VX(1,1:NIX+1)=VX(2,1:NIX+1) --> CLRFVX(1,1:NIX+1)=CLRFVX(2,1:NIX+1) --> Some mass flux along a-axis boundary --> top FREE boundary
        %VX(NIY+2,1:NIX+1)=VX(NIY+1,1:NIX+1) --> CLRFVX(NIY+2,1:NIX+1)=CLRFVX(NIY+1,1:NIX+1) --> Some mass flux along x-axis boundary --> bottom FREE boundary
        CLRFVX(j,i)=CLTM(j,i)*max(RFVX(j,i),0.0)+CLTM(j,i+1)*min(RFVX(j,i),0.0);%[kg/m^2/sec]
    end
end

%Recheck boundary condition explicitly
% for i=1:NIX+1
%     %Top NO SLIP boundary
%     CLRFVX(1,i)=-CLRFVX(2,i);
%     
%     %bottom NO SLIP boundary
%     CLRFVX(NIY+2,i)=-CLRFVX(NIY+1,i);
% end

for i=1:NIX+1
    %Top NO SLIP boundary
    CLRFVX(1,i)=-CLRFVX(2,i);
    
    %bottom FREE boundary
    CLRFVX(NIY+2,i)=CLRFVX(NIY+1,i);
end

%y-axis species liquid convection flux
CLRFVY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    %VY(1:NIY+1,1)=VY(1:NIY+1,2) --> CLRFVY(1:NIY+1,1)=CLRFVY(1:NIY+1,2) --> Some mass flux along y-axis boundary --> left FREE boundary
    %VY(1:NIY+1,NIX+2)=VY(1:NIY+1,NIX+1) --> CLRFVY(1:NIY+1,NIX+2)=CLRFVY(1:NIY+1,NIX+1) --> Some mass flux along y-axis boundary --> right FREE boundary
    
    %VY(1:NIY+1,1)=-VY(1:NIY+1,2) --> CLRFVY(1:NIY+1,1)=-CLRFVY(1:NIY+1,2) --> No mass flux along y-axis boundary --> left NO SLIP boundary
    %VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1) --> CLRFVY(1:NIY+1,NIX+2)=-CLRFVY(1:NIY+1,NIX+1) --> No mass flux along y-axis boundary --> right NO SLIP boundary
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> CLRFVY(1,1:NIX+2)=0.0 (top impermeable boundary)
        %VY(NIY+1,1:NIX+2)=0.0 --> CLRFVY(NIY+1,1:NIX+2)=0.0 (bottom impermeable boundary)
        CLRFVY(j,i)=CLTM(j,i)*max(RFVY(j,i),0.0)+CLTM(j+1,i)*min(RFVY(j,i),0.0);%[kg/m^2/sec]
    end
end

%Recheck boundary condition explicitly
for i=1:NIY+1
    CLRFVY(i,1)=CLRFVY(i,2);%left FREE boundary
    CLRFVY(i,NIX+2)=CLRFVY(i,NIX+1);%right FREE boundary
end

% for i=1:NIY+1
%     CLRFVY(i,1)=-CLRFVY(i,2);%left NO SLIP boundary
%     CLRFVY(i,NIX+2)=-CLRFVY(i,NIX+1);%right NO SLIP boundary
% end

%---------------------- Part One: CLRFVX CLRFVY ---------------------------

%------------------------ Part Two: Species flux --------------------------

%x-axis species liquid convection flux
CLQVX=zeros(NIY+2,NIX);%CL=CL, Q=flux, VX=VX
for i=1:NIX
    for j=1:NIY+2
        CLQVX(j,i)=-dtb*(CLRFVX(j,i+1)-CLRFVX(j,i))/dx(i);%[kg/m^3]
    end
    %CLQVX(1,:) and CLQVX(NIY+2,:) is useless, its value doesn't matter
end

%y-axis species liquid convection flux
CLQVY=zeros(NIY,NIX+2);%CL=CL, Q=flux, VY=VY
for i=1:NIX+2
    for j=1:NIY
        CLQVY(j,i)=-dtb*(CLRFVY(j+1,i)-CLRFVY(j,i))/dy(j);%[kg/m^3]
    end
    %CLQVY(:,1) and CLQVY(:,NIX+2) is useless, its value doesn't matter
end

%total species liquid convection flux
CLQVT=zeros(NIY,NIX);%CL=CL, Q=flux, V=velocity, T=total
CLQVT=CLQVX(2:NIY+1,1:NIX)+CLQVY(1:NIY,2:NIX+1);%[kg/m^3]

%------------------------ Part Two: Species flux --------------------------

%% =============== SOLID CONVECTION SPECIES FLUX =================


%% ============= RESIDUAL INTERNAL SPECIES STORAGE ===============

% KP=GeoConst.KSL;%solid/liquid partition coefficient both for major and trace elements
%NOTE: for major elements, distribution coefficients can be found Elkins-Tanton 2008 Supplementary; while for trace elements, partition coefficients
%can also be found in Elkins-Tanton 2008 Supplementary
% KP=zeros(NIY+2,NIX+2);%solid/liquid partition coefficient
% for i=1:NIX+2
%     for j=1:NIY+2
%         if(CLTM(j,i)<=CE)
%             KP(j,i)=0.12824+5.699124e-5*CLTM(j,i)*100.0+3.728777e-5*(CLTM(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
%         else
%             KP(j,i)=1.0;
%         end
%     end
% end

%ISS=zeros(NIY+2,NIX+2);%Residual Internal Species Storage

% for i=2:NIX+1
%     for j=2:NIY+1
%         %ISS(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
%         %NOTE: The above expression will give unphysical result when
%         %crystallization does not begin: CLTM does not change but RLTM will
%         %increase due to temperature drop!
%         
%         ISS.ANY(j,i)=ISS.ANY(j,i)+CLQVT.ANY(j-1,i-1)+CLQDT.ANY(j-1,i-1)+CSQVT.ANY(j-1,i-1)-CLFNS.ANY(j-1,i-1);%[kg/m^3]
%         %NOTE: the above expression is always valid, and used in MAGTFC.m
%     end
% end

end

