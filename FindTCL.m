function FindTCL(FLTM,RLTM,RSTM)
%This function calculates trace elements concentration  in liquid.
%Created on 2020-4-12

global NIX
global NIY
global dx
global dy
global dtb
global TCL
global DL
global Minors
global KSL
global dFS

phi=0.5;
FRCLX=zeros(NIY+2,NIX+2);
FRCLY=zeros(NIY+2,NIX+2);
CLQDX=zeros(NIY+2,NIX+1);%CL=CL, Q=flux, D=diffusion, X=x-axis
CLQDY=zeros(NIY+1,NIX+2);%CL=CL, Q=flux, D=diffusion, Y=y-axis
CLQDT=zeros(NIY,NIX);%CL=CL, Q=flux, D=diffusion, T=total
ISS=zeros(NIY+2,NIX+2);%internal species storage
FRLTM=zeros(NIY+2,NIX+2);%old time FL*RL
RSKPTM=zeros(NIY+2,NIX+2);%old time RS*KP
CLRFVX=zeros(NIY+2,NIX+1);
CLRFVY=zeros(NIY+1,NIX+2);
CLQVX=zeros(NIY+2,NIX);%CL=CL, Q=flux, VX=VX
CLQVY=zeros(NIY,NIX+2);%CL=CL, Q=flux, VY=VY
CLQVT=zeros(NIY,NIX);%CL=CL, Q=flux, V=velocity, T=total
m=length(fieldnames(TCL));
n=length(fieldnames(dFS));
dFST=0.0;

for k=1:m
%% .............========= DIFFUSION SPECIES FLUX =========.................

cmd=['TCL.',Minors{k}];
CLTM=eval(cmd);
%x-axis FL*RL*CL
% FRCLX=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        FRCLX(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
    end
end

%y-axis FL*RL*CL
% FRCLY=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        FRCLY(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
    end
end

%x-axis diffusion
% CLQDX=zeros(NIY+2,NIX+1);%CL=CL, Q=flux, D=diffusion, X=x-axis
for i=2:NIX
    %CL(1:NIY+2,1)=CL(1:NIY+2,2) --> CLQDX(1:NIY+2,1)=0.0 (left insulated boundary)
    %CL(1:NIY+2,NIX+1)=CL(1:NIY+2,NIX+2) --> CLQDX(1:NIY+2,NIX+1)=0.0 (right impermeable boundary)
    for j=1:NIY+2
        CLQDX(j,i)=2.0*dtb*DL*(FRCLX(j,i+1)-FRCLX(j,i))/((dx(i-1)+dx(i))*dx(i-1));%[kg/m^3]
    end
end

%y-axis diffusion
% CLQDY=zeros(NIY+1,NIX+2);%CL=CL, Q=flux, D=diffusion, Y=y-axis
for i=1:NIX+2
    for j=2:NIY
        %CL(1,1:NIX+2)=CL(2,1:NIX+2) --> CLQDY(1,1:NIX+2)=0.0 (top impermeable boundary)
        %CL(NIY+2,1:NIX+2)=CL(NIY+1,1:NIX+2) --> CLQDY(NIY+1,1:NIX+2)=0.0 (bottom impermeable boundary)
        CLQDY(j,i)=2.0*dtb*DL*(FRCLY(j+1,i)-FRCLY(j,i))/((dy(j-1)+dy(j))*dy(j-1));%[kg/m^3]
    end
end

%total species flux of each finite volume
% CLQDT=zeros(NIY,NIX);%CL=CL, Q=flux, D=diffusion, T=total
for i=1:NIX
    for j=1:NIY
        CLQDT(j,i)=CLQDX(j+1,i+1)-CLQDX(j+1,i)+CLQDY(j+1,i+1)-CLQDY(j,i+1);%[kg/m^3]
    end
end

%% .............========= INTERNAL SPECIES STORAGE =========...............

% ISS=zeros(NIY+2,NIX+2);%internal species storage
% FRLTM=zeros(NIY+2,NIX+2);%old time FL*RL
% RSKPTM=zeros(NIY+2,NIX+2);%old time RS*KP

for i=1:NIX+2
    for j=1:NIY+2
        FRLTM(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
        KSLValue=eval(['KSL.',Minors{m},'(j,i)']);
        RSKPTM(j,i)=RSTM(j,i)*KSLValue*CLTM(j,i);%[kg/m^3]
        dFS.OL(j,i)
        ISS(j,i)=FRLTM(j,i)-(1.0-phi)*dFS(j,i)*RSKPTM(j,i);%[kg/m^3]
    end
end

%% ...........========== CONVECTION SPECIES FLUX ===========...............

%---------------------- Part One: CLRFVX CLRFVY ---------------------------

%x-axis species convection flux
% CLRFVX=zeros(NIY+2,NIX+1);
for i=2:NIX
    %VX(1:NIY+2,1)=0.0 --> CLRFVX(1:NIY+2,1)=0.0 (left insulated boundary)
    %VX(1:NIY+2,NIX+1)=0.0 --> CLRFVX(1:NIY+2,NIX+1)=0.0 (right impermeable boundary)
    for j=1:NIY+2
        %VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> CLRFVX(1,1:NIX+1)=-CLRFVX(2,1:NIX+1) (top impermeable boundary)
        %VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1) --> CLRFVX(NIY+2,1:NIX+1)=-CLRFVX(NIY+1,1:NIX+1) (bottom impermeable boundary)
        CLRFVX(j,i)=CLTM(j,i)*max(RFVX(j,i),0.0)+CLTM(j,i+1)*min(RFVX(j,i),0.0);%[kg/m^2/sec]
    end
end

%y-axis species convection flux
% CLRFVY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    %VY(1:NIY+1,1)=FREE --> CLRFVY(1:NIY+1,1)==CLRFVY(1:NIY+1,2) (left slip boundary)
    %VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1) --> CLRFVY(1:NIY+1,NIX+2)=-CLRFVY(1:NIY+1,NIX+1) (right impermeable boundary)
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> CLRFVY(1,1:NIX+2)=0.0
        %VY(NIY+1,1:NIX+2)=0.0 --> CLRFVY(NIY+1,1:NIX+2)=0.0
        CLRFVY(j,i)=CLTM(j,i)*max(RFVY(j,i),0.0)+CLTM(j+1,i)*min(RFVY(j,i),0.0);%[kg/m^2/sec]
    end
end

for i=1:NIY+1
    CLRFVY(i,1)=CLRFVY(i,2);%left slip boundary
end

%---------------------- Part One: CLRFVX CLRFVY ---------------------------

%------------------------ Part Two: Species flux --------------------------

%x-axis species convection flux
% CLQVX=zeros(NIY+2,NIX);%CL=CL, Q=flux, VX=VX
for i=1:NIX
    for j=1:NIY+2
        CLQVX(j,i)=-dtb*(CLRFVX(j,i+1)-CLRFVX(j,i))/dx(i);%[kg/m^3]
    end
end

%y-axis species convection flux
% CLQVY=zeros(NIY,NIX+2);%CL=CL, Q=flux, VY=VY
for i=1:NIX+2
    for j=1:NIY
        CLQVY(j,i)=-dtb*(CLRFVY(j+1,i)-CLRFVY(j,i))/dy(j);%[kg/m^3]
    end
end

%total species convection flux
% CLQVT=zeros(NIY,NIX);%CL=CL, Q=flux, V=velocity, T=total
CLQVT=CLQVX(2:NIY+1,1:NIX)+CLQVY(1:NIY,2:NIX+1);%[kg/m^3]

%------------------------ Part Two: Species flux --------------------------
end

%New TCL
CLFN=zeros(NIY,NIX);%CL Factor of Next step
CLN=zeros(NIY+2,NIX+2);%CL of Next step
KPN=zeros(NIY+2,NIX+2);%KP of Next step

for i=1:NIX
    for j=NIY
                CLFN(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1)+phi*dFSTM(j+1,i+1)*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
                CLN(j+1,i+1)=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFN(j,i);%[kg/m^3 divided by kg/m^3 ==1]; ISS should update with NEW dFS ,but no NEW dFS available now, set old dFSTM==0.0 as new one
    end
end

end

