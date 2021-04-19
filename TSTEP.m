function  TimeInterval= TSTEP(FRCPSTM,FLTM,FSTM,RLTM,RSTM,CPLTM,KLTM,FKTM,MUOTM,VXTM,VYTM,step)%VSXTM,VSYTM
%To estimate time interval between n and n+1 step
%Created 2019-10-17
%Modified for crystal movement dt constraint 2020-6-30

%FRCPSTM [NIY+2,NIX+2]
%FLTM [NIY+2,NIX+2]
%RLTM [NIY+2,NIX+2]
%CPLTM [NIY+2,NIX+2]
%KLTM [NIY+2,NIX+2]
%FKTM [NIY+2,NIX+2]

%VXTM [NIY+2,NIX+1]
%VYTM [NIY+1,NIX+2]

global NIX
global NIY
global dx
global dy
global FSCR1

a=0.6;%to relax dt

%% ================= POSITIVE TEMPERATURE LIMIT ==================
IHS=zeros(NIY,NIX);%old T related Internal Heat Storage

%-------------------------------- TFO -------------------------------------
for i=1:NIX
    for j=1:NIY
        IHS(j,i)=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1);%[J/m^3/K]
    end
end

%-------------------------------- TJUG ------------------------------------
%x-axis velocity
%        +-------+
%        |       |
%    --> |   *   |  -->
%        |       |
%        +-------+
%NOTE: in fact, vx at loft and bottom faces should be estimated,
%but they can be interpolated; so only vx normal to vertical faces are
%calculated!

RFVX=zeros(NIY+2,NIX+1);%RL*FL*VX, VX faces
for i=2:NIX
    %VX(1:NIY+2,1)=0.0 --> RFVX(1:NIY+2,1)=0.0 (left impermeable boundary)
    %VX(1:NIY+2,NIX+1)=0.0 --> RFVX(1:NIY+2,NIX+1)=0.0 (right impermeable boundary)
    for j=2:NIY+1
        %VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> RFVX(1,1:NIX+1)=-RFVX(2,1:NIX+1) (top NO SLIP boundary)
        %VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1) --> RFVX(NIY+2,1:NIX+1)=-RFVX(NIY+1,1:NIX+1) (bottom NO SLIP boundary)
        
        %VX(1,1:NIX+1)=VX(2,1:NIX+1) --> RFVX(1,1:NIX+1)=RFVX(2,1:NIX+1) (top FREE boundary)
        %VX(NIY+2,1:NIX+1)=VX(NIY+1,1:NIX+1) --> RFVX(NIY+2,1:NIX+1)=RFVX(NIY+1,1:NIX+1) (bottom FREE boundary)
        
        RFVX(j,i)=VXTM(j,i)*0.5*(RLTM(j,i)+RLTM(j,i+1))*0.5*(FLTM(j,i)+FLTM(j,i+1));%[kg/m^2/sec]
    end
end

for i=1:NIX+1
    RFVX(1,i)=-RFVX(2,i);%top NO SLIP boundary
    RFVX(NIY+2,i)=RFVX(NIY+1,i);%bottom FREE boundary
end

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

RFVY=zeros(NIY+1,NIX+2);%RL*FL*VY, VY faces in ingot
for i=2:NIX+1
    %VY(1:NIY+1,1)=VY(1:NIY+1,2) --> RFVY(1:NIY+1,1)=RFVY(1:NIY+1,2) (left FREE boundary)
    %VY(1:NIY+1,NIX+2)=VY(1:NIY+1,NIX+1) --> RFVY(1:NIY+1,NIX+2)=RFVY(1:NIY+1,NIX+1) (right FREE boundary)
    
    %VY(1:NIY+1,1)=-VY(1:NIY+1,2) --> RFVY(1:NIY+1,1)=-RFVY(1:NIY+1,2) (left NO SLIP boundary)
    %VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1) --> RFVY(1:NIY+1,NIX+2)=-RFVY(1:NIY+1,NIX+1) (right NO SLIP boundary)
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> RFVY(1,1:NIX+2)=0.0 (top impermeable boundary)
        %VY(NIY+1,1:NIX+2)=0.0 --> RFVY(NIY+1,1:NIX+2)=0.0 (bottom impermeable boundary)
        RFVY(j,i)=VYTM(j,i)*0.5*(RLTM(j,i)+RLTM(j+1,i))*0.5*(FLTM(j,i)+FLTM(j+1,i));%[kg/m^2/sec]
    end
end

for i=1:NIY+1
    RFVY(i,1)=RFVY(i,2);%FREE left boundary
    RFVY(i,NIX+2)=RFVY(i,NIX+1);%FREE right boundary
end

CPLRFV=zeros(NIY,NIX);%CPL*RFVX, CPL*RFVY
for i=1:NIX
    for j=1:NIY
        CPLRFV(j,i)=CPLTM(j+1,i+1)*((max(RFVX(j+1,i+1),0.0)-min(RFVX(j+1,i),0.0))/dx(i)+(max(RFVY(j+1,i+1),0.0)-min(RFVY(j,i+1),0.0))/dy(j));%[J/K/m^3/sec]
    end
end

%------------- thermal conductivity on finite volume faces ----------------

%........................all faces, x-axis [W/m/K] ........................
KX=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    for j=1:NIY+2
        KX(j,i)=0.5*(FKTM(j,i)+FLTM(j,i)*KLTM(j,i)+FKTM(j,i+1)+FLTM(j,i+1)*KLTM(j,i+1));%[(FS*KS)+(FL*KL)]|(j+/-0.5,k); KX at T(1,1) and T(1,2) is not useful and so omitted
    end
end

%........................all faces, x-axis [W/m/K] ........................

%........................all faces, y-axis [W/m/K] ........................
KY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        KY(j,i)=0.5*(FKTM(j,i)+FLTM(j,i)*KLTM(j,i)+FKTM(j+1,i)+FLTM(j+1,i)*KLTM(j+1,i));%[(FS*KS)+(FL*KL)]|(j,k+/-0.5); KY at T(1,1) and T(2,1) is not useful and so omitted
    end
end

%........................all faces, y-axis [W/m/K] ........................

%------------- thermal conductivity on finite volume faces ----------------

FKD=zeros(NIY,NIX);%(FS*KS+FL*KL)/(dx or dy)
for i=2:NIX-1
    for j=2:NIY-1
        FKD(j,i)=2.0*KX(j+1,i+1)/((dx(i)+dx(i+1))*dx(i))+2.0*KX(j+1,i)/((dx(i)+dx(i-1))*dx(i))+2.0*KY(j+1,i+1)/((dy(j)+dy(j+1))*dy(j))+2.0*KY(j,i+1)/((dy(j)+dy(j-1))*dy(j));%[W/m^3/K]
    end
end

for i=2:NIX-1
    FKD(1,i)=2.0*KX(2,i)/((dx(i-1)+dx(i))*dx(i))+2.0*KX(2,i+1)/((dx(i)+dx(i+1))*dx(i))+2.0*KY(1,i+1)/((dy(1)+dy(1))*dy(1))+2.0*KY(2,i+1)/((dy(1)+dy(2))*dy(1));%top boundary [W/m^3/K]
    FKD(NIY,i)=2.0*KX(NIY+1,i)/((dx(i-1)+dx(i))*dx(i))+2.0*KX(NIY+1,i+1)/((dx(i)+dx(i+1))*dx(i))+2.0*KY(NIY,i+1)/((dy(NIY-1)+dy(NIY))*dy(NIY))+2.0*KY(NIY+1,i+1)/((dy(NIY)+dy(NIY))*dy(NIY));%bottom boundary [W/m^3/K]
end

for i=2:NIY-1
    FKD(i,1)=2.0*KX(i+1,1)/((dx(1)+dx(1))*dx(1))+2.0*KX(i+1,2)/((dx(1)+dx(2))*dx(1))+2.0*KY(i,2)/((dy(i-1)+dy(i))*dy(i))+2.0*KY(i+1,2)/((dy(i)+dy(i+1))*dy(i));%left boundary [W/m^3/K]
    FKD(i,NIX)=2.0*KX(i+1,NIX)/((dx(NIX-1)+dx(NIX))*dx(NIX))+2.0*KX(i+1,NIX+1)/((dx(NIX)+dx(NIX))*dx(NIX))+2.0*KY(i,NIX+1)/((dy(i-1)+dy(i))*dy(i))+2.0*KY(i+1,NIX+1)/((dy(i)+dy(i+1))*dy(i));%right boundary [W/m^3/K]
end

FKD(1,1)=2.0*KX(2,1)/((dx(1)+dx(1))*dx(1))+2.0*KX(2,2)/((dx(2)+dx(1))*dx(1))+2.0*KY(1,2)/((dy(1)+dy(1))*dy(1))+2.0*KY(2,2)/((dy(1)+dy(2))*dy(1));%top-left point
FKD(NIY,1)=2.0*KX(NIY+1,1)/((dx(1)+dx(1))*dx(1))+2.0*KX(NIY+1,2)/((dx(2)+dx(1))*dx(1))+2.0*KY(NIY,2)/((dy(NIY-1)+dy(NIY))*dy(NIY))+2.0*KY(NIY+1,2)/((dy(NIY)+dy(NIY))*dy(NIY));%bottom-left point
FKD(1,NIX)=2.0*KX(2,NIX)/((dx(NIX-1)+dx(NIX))*dx(NIX))+2.0*KX(2,NIX+1)/((dx(NIX)+dx(NIX))*dx(NIX))+2.0*KY(1,NIX+1)/((dy(1)+dy(1))*dy(1))+2.0*KY(2,NIX+1)/((dy(2)+dy(1))*dy(1));%right-top point
FKD(NIY,NIX)=2.0*KX(NIY+1,NIX)/((dx(NIX-1)+dx(NIX))*dx(NIX))+2.0*KX(NIY+1,NIX+1)/((dx(NIX)+dx(NIX))*dx(NIX))+2.0*KY(NIY,NIX+1)/((dy(NIY-1)+dy(NIY))*dy(NIY))+2.0*KY(NIY+1,NIX+1)/((dy(NIY)+dy(NIY))*dy(NIY));%right-bottom point

dtD=zeros(NIY,NIX);%dt=dt, D=diffusion
for i=1:NIX
    for j=1:NIY
        dtD(j,i)=IHS(j,i)/(CPLRFV(j,i)+FKD(j,i));%IHS=[NIY,NIX] [sec]
    end
end

dt1=min(min(dtD));%diffusion time limit

%% ====================== DIFFUSION LIMIT ========================
%Courant£­Friedrichs£­Lewy limit: dt<=dx^2/(3*kap)
Lmin=min(min(dx),min(dy));%minimum length scale [m]
Kmax=max(max(KX(:)),max(KY(:)));%maximum thermal conductivity [W/m/K]
RCPmin=min(min(FRCPSTM+FLTM.*RLTM.*CPLTM));%minimum rho*CPS [J/m^3/K]
dt2=Lmin^2/(3.0*Kmax/RCPmin);

%% ===================== CONVECTION LIMIT =========================
%liquid x-axis velocity limit
VXL=zeros(NIY+2,NIX+1);%L==limit
for i=2:NIX
    %VX(1:NIY+2,1)=0.0 --> dt -> Inf
    %VX(1:NIY+2,NIX+1)=0.0 --> dt -> Inf
    for j=1:NIY+2
        %VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> dt < Inf (top NO SLIP boundary)
        %VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1)--> dt < Inf (bottom NO SLIP boundary)
        
        %VX(1,1:NIX+1)=VX(2,1:NIX+1) --> dt < Inf (top FREE boundary)
        %VX(NIY+2,1:NIX+1)=VX(NIY+1,1:NIX+1) --> dt < Inf (bottom FREE boundary)
        VXL(j,i)=(dx(i-1)+dx(i))/(2.0*abs(VXTM(j,i)));%[sec]
    end
end
VXL(1:NIY+2,1)=dx(1)./(abs(VXTM(1:NIY+2,1)));
VXL(1:NIY+2,NIX+1)=dx(NIX)./abs(VXTM(1:NIY+2,NIX+1));
VXMIN=min(min(VXL(2:NIY+1,2:NIX)));

%liquid y-axis velocity limit
VYL=zeros(NIY+1,NIX+2);%L==limit
for i=1:NIX+2
    %VY(1:NIY+1,1)=VY(1:NIY+1,2) --> dt < Inf (left FREE boundary)
    %VY(1:NIY+1,NIX+2)=VY(1:NIY+1,NIX+1) --> dt < Inf (right FREE boundary)
    
    %VY(1:NIY+1,1)=-VY(1:NIY+1,2) --> dt < Inf (left NO SLIP boundary)
    %VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1) --> dt < Inf (right NO SLIP boundary)
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> dt -> Inf (top impermeable boundary)
        %VY(NIX+1,1:NIX+2)=0.0 --> dt -> Inf (bottom impermeable boundary)
        VYL(j,i)=(dy(j-1)+dy(j))/(2.0*abs(VYTM(j,i)));%[sec]
    end
end
VYL(1,1:NIX+2)=dy(1)./abs(VYTM(1,1:NIX+2));
VYL(NIY+1,1:NIX+2)=dy(NIY)./abs(VYTM(NIY+1,1:NIX+2));
VYMIN=min(min(VYL(2:NIY,2:NIX+1)));

dt3=min(VXMIN,VYMIN);%convection time limit

%% ======================== TIME LIMIT ===========================
if(step<=20)
    TimeInterval=a*dt1/100.0;
end
if(20<step&&step<=30)
TimeInterval=a*dt1/200.0;

end
if(30<step&&step<=150)
    TimeInterval=a*dt1/100.0;
end

if(step>150)
    TimeInterval=a*min([dt1,dt2,dt3]);
end

%TimeInterval=min([dt1,dt2,dt3]);
%10 m x 10 m cavity with 20 x 20 grid has maximum 10^3 sec time interval, so set as standard

end

