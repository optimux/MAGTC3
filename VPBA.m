function [VX,VY,P,RESM] = VPBA(VXOTM,VYOTM,PDTM,MUOTM,MUNTM,FLOTM,FLNTM,RLOTM,RLNTM,RSOTM,RSNTM,dFSTM,RLOGTM,RLNGTM)
%To solve V-P momentum conservation
%Created on 2019-10-20

%Modified for iteration logical 2019-11-10

%Modified for boundary condition settings 2019-11-13

%Modified for old terms corrected from new terms 2020-2-14

%Modified for mass check 2020-2-16

%Modified for dynamic pressure and gravity term 2020-3-26

%Modified for crystal movement 2020-7-9

%Add full comments on boundary condition 2020-8-12

%O==Old
%N==New
%BA-->Boussinesq Approximation

global NIX
global NIY
global dx
global dy
global g
global RL0
global dtb
global P0
global FSCR1
global FS
global RS
global FSO
global RSO
global dFS
global K0

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%compound variables. These packages can be utilized repeteadly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================
%We add the following boundary cells. Dummy coefficients FLNTM(1,:)==0.0 and FLOTM(1,:)==0.0 because real top boundary is cold and pure solid. FLNTM(1,:)=0.0
%and FLOTM(1,:) have been used in pgysical properties calculation, so we modified them here temporarily to make a ghost cell to make boundary conditions
%resonable.
FLNTM(1,:)=FLNTM(2,:);
FLOTM(1,:)=FLOTM(2,:);

FLErr=1.0e-16;%shreshold of error for determining FL=0.0, i.e., if FLNTM(j,i)<=FLErr, then FLNTM(j,i)=0.0. In most cases, FLO(j,i)<1.0e-4 then FLN==0.0!
VErr=1.0e-16;%shreshold of error for determining VX, VY=0.0

%% ======================= DFLP COEFFICIENT MATRIX =========================

%--------------------------- Permeability K -------------------------------
FSNTM=1.0-FLNTM;
FSOTM=1.0-FLOTM;
%Permeability in Porus Region and Free-Moving Region [m^2]
KFLTM=zeros(NIY+2,NIX+2);
% Xudaming 1991
% for i=1:NIX+2
%     for j=1:NIY+2
%         if(FLNTM(j,i)>=1.0/3.0)
%             KFLTM(j,i)=2.6e-5*(1.923e-2*FLNTM(j,i)^2+(4.0+3.0*FSNTM(j,i)-3.0*real(sqrt(FSNTM(j,i)*(8.0-3.0*FSNTM(j,i)))))/FSNTM(j,i));%[mm^2]
%         else
%             KFLTM(j,i)=5.0e-7*FLNTM(j,i)^2;%[mm^2]
%         end
%     end
% end
% KFLTM=KFLTM*10^-6;%[m^2]

% Blake-Kozeny-Carman
for i=1:NIX+2
    for j=1:NIY+2
        Ga=1.0;%(0.5+atan(100.0*(1.0-FSCR1-FLNTM(j,i)))/pi)^-5;%rheology transition factor is 5, Ga should be omitted here since no solid motion.
        KFLTM(j,i)=Ga*K0*(FLNTM(j,i)^3/(1.0-FLNTM(j,i))^2);
    end
end

%--------------------------- Permeability K -------------------------------

%------------------- mu*FL/(RL*K) in denomenator --------------------------
UFRKY=zeros(NIY+1,NIX+2);%Y==y-axis
for i=1:NIX+2
    for j=1:NIY+1
        %IMPORTANT NOTE: at first, 1.0e-16 is added to prevent devision by 0, this
        %helps when FL=0.0; but this addition is discarded and emplaced by
        %the following.
        UFRKY(j,i)=0.5*(MUNTM(j,i)+MUNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i))/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(KFLTM(j,i)+KFLTM(j+1,i)));%[1/sec]
        
        %solid cell at (j,i), set UFRKY as 0.0 mandatorily
        if((FLNTM(j,i)<=FLErr)&&(FLNTM(j+1,i)<=FLErr))
            UFRKY(j,i)=0.0;
        end
        
        %For free-moving region, no Darcy damping force, set UFRKY as 0.0 mandatorily
        if(0.5*(FSNTM(j,i)+FSNTM(j+1,i))<=FSCR1)
            UFRKY(j,i)=0.0;
        end
        
    end
end

UFRKX=zeros(NIY+2,NIX+1);%X==x-axis
for i=1:NIX+1
    for j=1:NIY+2
        %IMPORTANT NOTE: at first, 1.0e-16 is added to prevent devision of 0, this
        %helps when FL=0.0; but this addition is discarded and emplaced by
        %the following.
        UFRKX(j,i)=0.5*(MUNTM(j,i)+MUNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1))/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(KFLTM(j,i)+KFLTM(j,i+1)));%[1/sec]
        
        %solid cell at (j,i), set UFRKX as 0.0 mandatorily
        if((FLNTM(j,i)<=FLErr)&&(FLNTM(j,i+1)<=FLErr))
            UFRKX(j,i)=0.0;
        end
        
        %For free-moving region, no Darcy damping force, set UFRKX as 0.0 mandatorily
        if(0.5*(FSNTM(j,i)+FSNTM(j,i+1))<=FSCR1)
            UFRKX(j,i)=0.0;
        end
        
    end
end

%------------------- mu*FL/(RL*K) in denomenator --------------------------
%x-axis parameter in matrix
AX=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        AX(j,i)=2.0*dy(j)/((dx(i-1)+dx(i))*(1.0+dtb*UFRKX(j+1,i)));%main domain [1]
    end
end

for i=1:NIY
    AX(i,1)=2.0*dy(i)/(2.0*dx(1)*(1.0+dtb*UFRKX(i+1,1)));%1st column [1]
    AX(i,NIX+1)=2.0*dy(i)/(2.0*dx(NIX)*(1.0+dtb*UFRKX(i+1,NIX+1)));%last column [1]
end

%y-axis parameter in matrix
AY=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        AY(j,i)=2.0*dx(i)/((dy(j)+dy(j-1))*(1.0+dtb*UFRKY(j,i+1)));%main domain [1]
    end
end

for i=1:NIX
    AY(1,i)=2.0*dx(i)/(2.0*dy(1)*(1.0+dtb*UFRKY(1,i+1)));%1st row [1]
    AY(NIY+1,i)=2.0*dx(i)/(2.0*dy(NIY)*(1.0+dtb*UFRKY(NIY+1,i+1)));%last row [1]
end

%=========================== [A] MATRIX ===============================
%NOTE: In [A][X]=[B], [A], [X] and [B] have dimension of [NIY]*[NIX]; to solve this 2D
%matrix, we reshape [X] and [B] into vector of length [NIX*NIY], this will give a
%large sparse matrix of [A] of dimension [NIX*NIY]*[NIX*NIY]

k=0;
%coefficient matrix, a(j,k), a(j+0.5,k), a(j-0.5,k), a(j,k-0.5), a(j,k+0.5)
A=zeros(NIY*NIX,NIX*NIY);
for i=2:NIX-1
    for j=2:NIY-1
        k=(i-1)*NIY+j;
        A(k,k-NIY)=-AX(j,i);%[1]
        A(k,k-1)=-AY(j,i);%[1]
        A(k,k)=AX(j,i)+AY(j,i)+AY(j+1,i)+AX(j,i+1);%[1]
        A(k,k+1)=-AY(j+1,i);%[1]
        A(k,k+NIY)=-AX(j,i+1);%[1]
    end
end

%special boundary condition at [1,1]
% A(1,1)=AX(1,2)+AY(2,1);%[1]
% A(1,2)=-AY(2,1)*0.0;%[1]
% A(1,NIY+1)=-AX(1,2)*0.0;%[1]
% A(1,1)=1.0;

A(1,1)=AX(1,2)+AY(2,1);%[1]
A(1,2)=-AY(2,1);%[1]
A(1,NIY+1)=-AX(1,2);%[1]

%special boundary condition at [NIY,1]
A(NIY,NIY)=AX(NIY,2)+AY(NIY,1);%[1]
A(NIY,NIY-1)=-AY(NIY,1);%[1]
A(NIY,2*NIY)=-AX(NIY,2);%[1]

%special boundary condition at [1,NIX]
A((NIX-1)*NIY+1,(NIX-1)*NIY+1)=AX(1,NIX)+AY(2,NIX);%[1]
A((NIX-1)*NIY+1,(NIX-1)*NIY+2)=-AY(2,NIX);%[1]
A((NIX-1)*NIY+1,(NIX-2)*NIY+1)=-AX(1,NIX);%[1]

%special boundary condition at [NIY,NIX]
A(NIX*NIY,NIX*NIY)=AX(NIY,NIX)+AY(NIY,NIX);%[1]
A(NIX*NIY,NIX*NIY-1)=-AY(NIY,NIX)*0.0;%[1]
A(NIX*NIY,(NIX-1)*NIY)=-AX(NIY,NIX)*0.0;%[1]
A(NIX*NIY,NIX*NIY)=1.0;

%Left boundary
for j=2:NIY-1
    A(j,j)=AY(j,1)+AY(j+1,1)+AX(j,2);%[1]
    A(j,j-1)=-AY(j,1);%[1]
    A(j,j+1)=-AY(j+1,1);%[1]
    A(j,j+NIY)=-AX(j,2);%[1]
end

%Right boundary
for i=2:NIY-1
    j=(NIX-1)*NIY+i;
    A(j,j)=AY(i,NIX)+AY(i+1,NIX)+AX(i,NIX);%[1]
    A(j,j-1)=-AY(i,NIX);%[1]
    A(j,j+1)=-AY(i+1,NIX);%[1]
    A(j,j-NIY)=-AX(i,NIX);%[1]
end

%Top boundary
for i=2:NIX-1
    j=(i-1)*NIY+1;
    A(j,j)=AX(1,i)+AY(2,i)+AX(1,i+1);%[1]
    A(j,j-NIY)=-AX(1,i);%[1]
    A(j,j+1)=-AY(2,i);%[1]
    A(j,j+NIY)=-AX(1,i+1);%[1]
end

%Bottom boundary
for i=2:NIX-1
    j=i*NIY;
    A(j,j)=AX(NIY,i)+AY(NIY,i)+AX(NIY,i+1);%[1]
    A(j,j-NIY)=-AX(NIY,i);%[1]
    A(j,j-1)=-AY(NIY,i);%[1]
    A(j,j+NIY)=-AX(NIY,i+1);%[1]
end
%=========================== [A] MATRIX ===============================

%% ===================== RFVXTM RFVYTM =======================
%NOTE: DFLP at 4 physical walls or boundaries are useless since coefficients are set zero there, thus DFLP has onlt NIX*NIY elements.
DFLP=zeros(NIY*NIX,1);%consider delta(FL*PD) or delta(PD) as one variable which has (NIX*NIY) elements in one column [X]
B=zeros(NIY*NIX,1);%set as [B] which has NIX*NIY elements in one column so that [A][X]=[B]
DRFVX=zeros(NIY+2,NIX+1);
DRFVY=zeros(NIY+1,NIX+2);
VX=zeros(NIY+2,NIX+1);%returned liquid x-axis velocity
VY=zeros(NIY+1,NIX+2);%returned y-axis liquid velocity

%Old step Effective (bulk) dynamic viscosity, I==integer [Pa.sec]
MUOEI=zeros(NIY+2,NIX+2);
MUOEI=VisMix(FSOTM,MUOTM);
% for i=1:NIX+2
%     for j=1:NIY+2
%         MUOEI(j,i)=VisMix(FSOTM,MUOTM);
%         %if(FSOTM<FSCR1)
%         %MUOEI(j,i)=f(FSOTM,MUOTM);
%         %end
%     end
% end

%Old step Effective (bulk) dynamic viscosity, H==half grid [Pa.sec]
MUOEH=zeros(NIY+1,NIX+1);
FSQM=zeros(NIY+1,NIX+1);%FS Quarter Mean
MUQM=zeros(NIY+1,NIX+1);%MU Quarter Mean
for i=1:NIX+1
    for j=1:NIY+1
        FSQM(j,i)=0.25*(FSOTM(j,i)+FSOTM(j+1,i)+FSOTM(j,i+1)+FSOTM(j+1,i+1));
        MUQM(j,i)=0.25*(MUOTM(j,i)+MUOTM(j+1,i)+MUOTM(j,i+1)+MUOTM(j+1,i+1));
        %         if(F<FSCR1)
        %             MUOEH(j,i)=VisMix(FSOTM,MUOTM);
        %         end
    end
end
MUOEH=VisMix(FSQM,MUQM);

%--------------------------- 1. Old RFVX ---------------------------------
%RFVX(j+0.5,k)
RFVX=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    for j=1:NIY+2
        RFVX(j,i)=0.5*(RLOTM(j,i)+RLOTM(j,i+1))*0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i);%[kg/m^2/sec]
    end
end
%RFVX(1:NIY+2,1)=0.0 --> Left impermeable boundary
%RFVX(1:NIY+2,NIX+1)=0.0 --> Right impermeable boundary

%RFVX(1,1:NIX+1)=-RFVX(2,1:NIX+1) --> Top NO SLIP boundary
%RFVX(NIY+2,1:NIX+1)=-RFVX(NIY+1,1:NIX+1) --> Bottom NO SLIP boundary

%RFVX(1,1:NIX+1)=RFVX(2,1:NIX+1) --> Top FREE boundary
%RFVX(NIY+2,1:NIX+1)=RFVX(NIY+1,1:NIX+1) --> Bottom FREE boundary
%--------------------------- 1. Old RFVX ---------------------------------

%--------------------------- 2. Old RFVY ---------------------------------
%RFVY(j,k+0.5)
RFVY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        RFVY(j,i)=0.5*(RLOTM(j,i)+RLOTM(j+1,i))*0.5*(FLOTM(j,i)+FLOTM(j+1,i))*VYOTM(j,i);%[kg/m^2/sec]
    end
end
%RFVY(1,1:NIX+2)=0.0 --> Top impermeable boundary
%RFVY(NIY+1,1:NIX+2)=0.0 --> Bottom impermeable boundary

%RFVY(1:NIY+1,1)=RFVY(1:NIY+1,2) --> Left FREE boundary
%RFVY(1:NIY+1,NIX+2)=RFVY(1:NIY+1,NIX+1) --> Right FREE boundary

%RFVY(1:NIY+1,1)=-RFVY(1:NIY+1,2) --> Left NO SLIP boundary
%RFVY(1:NIY+1,NIX+2)=-RFVY(1:NIY+1,NIX+1) --> Right NO SLIP boundary
%--------------------------- 2. Old RFVY ---------------------------------

%----------------------- 3. RFVSX Increment ------------------------------

%----------------------- 3. RFVSX Increment ------------------------------

%----------------------- 4. RFVSY Increment ------------------------------

%----------------------- 4. RFVSY Increment ------------------------------

%%==================== VELOCITY APPROXIMATION PARTS ======================
%------------------------- 5. RFVX in upwind  -------------------------------
%RFVX(j,k)=RFVXI
RFVXI=zeros(NIY+2,NIX+2);%I==integer grid (main grid)
for i=2:NIX+1
    for j=2:NIY+1
        RFVXI(j,i)=RLOTM(j,i)*FLOTM(j,i)*0.5*(VXOTM(j,i-1)+VXOTM(j,i));%[kg/m^2/sec]
    end
end
for i=2:NIY+1
    %VX(1:NIY+2,1)=0.0 --> Left impermeable boundary
    RFVXI(i,1)=-RFVXI(i,2);
    
    %VX(1:NIY+2,NIX+1)=0.0 --> Right impermeable boundary
    RFVXI(i,NIX+2)=-RFVXI(i,NIX+1);
end

%Top & Bottom NO SLIP boundary
for i=2:NIX+1
    %VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> No mass flow along x-axis --> top NO SLIP boundary
    RFVXI(1,i)=-RFVXI(2,i);
    
    %VX(NIY+2,1:NIX+1)=VX(NIY+1,1:NIX+1) --> Some mass flow along x-axis --> bottom FREE boundary
    RFVXI(NIY+2,i)=RFVXI(NIY+1,i);
end

% %Top & Bottom FREE boundary
% for i=2:NIX+1
%     %VX(1,1:NIX+1)=VX(2,1:NIX+1) --> Some mass flow along x-axis --> top FREE boundary
%     RFVXI(1,i)=RFVXI(2,i);
%     
%     %VX(NIY+2,1:NIX+1)=VX(NIY+1,1:NIX+1) --> Some mass flow along x-axis --> bottom FREE boundary
%     RFVXI(NIY+2,i)=RFVXI(NIY+1,i);
% end

%RFVX(j+0.5,k+0.5)=RFVXH
RFVXH=zeros(NIY+1,NIX+1);%H==half grid (cross grid)
for i=1:NIX+1
    for j=1:NIY+1
        RFVXH(j,i)=0.25*(RLOTM(j,i)+RLOTM(j+1,i)+RLOTM(j,i+1)+RLOTM(j+1,i+1))*0.25*(FLOTM(j,i)+FLOTM(j+1,i)+FLOTM(j,i+1)+FLOTM(j+1,i+1))*0.5*(VXOTM(j+1,i)+VXOTM(j,i));%[kg/m^2/sec]
    end
end
%VX(1:NIY+2,1)=0.0 --> RFVXH(1:NIY+1,1)=0.0 --> Left impermeable boundary
%VX(1:NIY+2,NIX+1)=0.0 --> RFVXH(1:NIY+1,NIX+1)=0.0 --> Right impermeable boundary

%VX(1,1:NIX+1)=VX(2,1:NIX+1) --> RFVXH(1,1:NIX+1)!=0.0 --> Top FREE boundary
%VX(NIY+2,1:NIX+1)=VX(NIX+1,1:NIX+1) --> RFVXH(NIY+1,1:NIX+1)!=0.0 --> Bottom FREE boundary

%VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> RFVXH(1,1:NIX+1)=0.0 --> Top NO SLIP boundary
%VX(NIY+2,1:NIX+1)=-VX(NIX+1,1:NIX+1) --> RFVXH(NIY+1,1:NIX+1)=0.0 --> Bottom NO SLIP boundary

%------------------------- 5. RFVX in upwind --------------------------------

%------------------------- 6. RFVSX in upwind -------------------------------

%------------------------- 6. RFVSX in upwind -------------------------------

%------------------------- 7. RFVY in upwind --------------------------------
%RFVY(j,k)=RFVYI
RFVYI=zeros(NIY+2,NIX+2);%I==integer grid (main grid)
for i=2:NIX+1
    for j=2:NIY+1
        RFVYI(j,i)=RLOTM(j,i)*FLOTM(j,i)*0.5*(VYOTM(j-1,i)+VYOTM(j,i));%[kg/m^2/sec]
    end
end

for i=2:NIX+1
    %VY(1,1:NIX+2)=0.0 --> Top impermeable boundary
    RFVYI(1,i)=-RFVYI(2,i);
    
    %VY(NIY+1,1:NIX+2)=0.0 --> Bottom impermeable boundary
    RFVYI(NIY+2,i)=-RFVYI(NIY+1,i);
end

%Left & Right NO SLIP boundary
% for i=2:NIY+1
%     %VY(1:NIY+1,1)=-VY(1:NIY+1,2) --> left NO SLIP boundary
%     RFVYI(i,1)=-RFVYI(i,2);
%
%     %VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1) --> right NO SLIP boundary
%     RFVYI(i,NIX+2)=-RFVYI(i,NIX+1);
% end

%Left & Right FREE boundary
for i=2:NIY+1
    %VY(1:NIY+1,1)=VY(1:NIY+1,2) --> left FREE boundary
    RFVYI(i,1)=RFVYI(i,2);
    
    %VY(1:NIY+1,NIX+2)=VY(1:NIY+1,NIX+1) --> right FREE boundary
    RFVYI(i,NIX+2)=RFVYI(i,NIX+1);
end

%RFVY(j+0.5,k+0.5)=RFVYH
RFVYH=zeros(NIY+1,NIX+1);%H==half grid (cross grid)
for i=1:NIX+1
    for j=1:NIY+1
        RFVYH(j,i)=0.25*(RLOTM(j,i)+RLOTM(j+1,i)+RLOTM(j,i+1)+RLOTM(j+1,i+1))*0.25*(FLOTM(j,i)+FLOTM(j+1,i)+FLOTM(j,i+1)+FLOTM(j+1,i+1))*0.5*(VYOTM(j,i+1)+VYOTM(j,i));%[kg/m^2/sec]
    end
end
%VY(1,1:NIX+2)=0.0 --> RFVYH(1,1:NIX+1)=0.0 --> Top impermeable boundary
%VY(NIY+1,1:NIX+2)=0.0 --> RFVYH(NIY+1,1:NIX+1)=0.0 --> Bottom impermeable boundary

%VY(1:NIY+1,1)=VY(1:NIY+1,2) --> RFVYH(1;NIY+1,1)!=0.0 --> Left FREE boundary
%VY(1:NIY+1,NIX+2)=VY(1:NIY+1,NIX+1) --> RFVYH(1;NIY+1,NIX+1)!=0.0 --> Right FREE boundary 

%VY(1:NIY+1,1)=-VY(1:NIY+1,2) --> RFVYH(1;NIY+1,1)=0.0 --> Left NO SLIP boundary
%VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1) --> RFVYH(1;NIY+1,NIX+1)=0.0 --> Right NO SLIP boundary 

%------------------------- 7. RFVY in upwind --------------------------------

%------------------------- 8. RFVSY in upwind -------------------------------

%------------------------- 8. RFVSY in upwind -------------------------------

%----------------------- 9. liquid x-momentum -------------------------------
%Upwind scheme for liquid x-axis momentum
VVXI=zeros(NIY+2,NIX+2);%I==integer grid (main grid)
for i=2:NIX+1
    for j=1:NIY+2
        VVXI(j,i)=VXOTM(j,i-1)*max(RFVXI(j,i),0.0)+VXOTM(j,i)*min(RFVXI(j,i),0.0);%including VVX1, VVX2 [kg/m/sec^2]
    end
end

for i=2:NIY+1
    %RFVXI(i,1)=RLOTM(i,1)*FLOTM(i,1)*0.5*(VXOTM(i,1)+(-VXOTM(i,2))), we set VXOTM(1:NIY+2,0)=-VXOTM(1:NIY+2,2) so that VXOTM(1:NIY+2,1)=0,
    %RFVXI(i,1)=-RFVXI(i,2) --> VVXI(2:NIY+1,1)=VVXI(2:NIY+1,2) --> VXCONV(i-1,1)=0.0 --> Left impermeable boundary 
    VVXI(i,1)=-VXOTM(i,2)*max(-RLOTM(i,1)*FLOTM(i,1)*0.5*VXOTM(i,2),0.0)+0.0;
    
    %RFVXI(i,NIX+2)=RLOTM(i,NIX+2)*FLOTM(i,NIX+2)*0.5*(VXOTM(i,NIX+1)+(-VXOTM(i,NIX))), we set VXOTM(1:NIY+2,NIX+2)=-VXOTM(1:NIY+2,NIX) so that
    %VXOTM(i,NIX+1)=0.0, RFVXI(i,NIX+2)=-RFVXI(i,NIX+1) --> VVXI(2:NIY+1,NIX+2)=VVXI(2:NIY+1,NIX+1) --> VXCONV(i-1,NIX+1)=0.0 --> Right impermeable boundary 
    VVXI(i,NIX+2)=0.0-VXOTM(i,NIX)*min(-RLOTM(i,NIX+2)*FLOTM(i,NIX+2)*0.5*VXOTM(i,NIX),0.0);
end

%Upwind scheme for y-axis liquid momentum
VVXH=zeros(NIY+1,NIX+1);%H==half grid (cross grid)
for i=1:NIX+1
    for j=1:NIY+1
        VVXH(j,i)=VXOTM(j,i)*max(RFVYH(j,i),0.0)+VXOTM(j+1,i)*min(RFVYH(j,i),0.0);%[kg/m/sec^2]
    end
end

%VY(1,1:NIX+2)=0.0 --> RFVYH(1,1:NIX+1)=0.0 --> VVXH(1,1:NIX+1)=0.0 --> Top impermeable boundary
%VY(NIY+1,1:NIX+2)=0.0 --> RFVYH(NIY+1,1:NIX+1)=0.0 --> VVXH(NIY+1,1:NIX+1)=0.0 --> Bottom impermeable boundary

%VX(1:NIY+2,1)=0.0 --> VVXH(1:NIY+1,1)=0.0 --> Left impermeable boundary
%VX(1:NIY+2,NIX+1)=0.0 --> VVXH(1;NIY+1,NIX+1)=0.0 --> Right impermeable boundary 

%----------------------- 9. liquid x-momentum -------------------------------

%----------------------- 10. solid x-momentum -------------------------------

%----------------------- 10. solid x-momentum -------------------------------

%----------------------- 11. liquid y-flow momentum -------------------------
%Upwind scheme for liquid y-axis momentum
VVYI=zeros(NIY+2,NIX+2);%I==integer grid (main grid)
for i=1:NIX+2
    for j=2:NIY+1
        VVYI(j,i)=VYOTM(j-1,i)*max(RFVYI(j,i),0.0)+VYOTM(j,i)*min(RFVYI(j,i),0.0);%including VVY3, VVY4 [kg/m/sec^2]
    end
end

for i=2:NIX+1
    %RFVYI(1,i)=RLOTM(1,i)*FLOTM(1,i)*0.5*(VYOTM(0,i)+VYOTM(1,i)), we set VYOTM(0,1:NIX+2)=-VYOTM(2,1:NIX+2) so that VYOTM(1,1:NIX+2)=0.0,
    %RFVYI(1,i)=-RFVYI(2,i) --> VVYI(1,2:NIX+1)=VVYI(2,2:NIX+1) --> VYCONV(1,i-1)=0.0 --> Top impermeable boundary
    VVYI(1,i)=-VYOTM(2,i)*max(-RLOTM(1,i)*FLOTM(1,i)*0.5*VYOTM(2,i),0.0)+0.0;
    
    %RFVYI(NIY+2,i)=RLOTM(NIY+2,i)*FLOTM(NIY+2,i)*0.5*(VYOTM(NIY+2,i)+VYOTM(NIY+1,i)), we set VYOTM(NIY+2,1:NIX+2)=-VYOTM(NIY,1:NIX+2) so that
    %VYOTM(NIY+1,1:NIX+2)=0.0, RFVYI(NIY+2,i)=-RFVYI(NIY+1,i) --> VVYI(NIY+2,2:NIX+1)=VVYI(NIY+1,2:NIX+1) --> VYCONV(NIY+1,i-1)=0.0 --> Bottom impermeable
    %boundary
    VVYI(NIY+2,i)=0.0-VYOTM(NIY,i)*min(-RLOTM(NIY+2,i)*FLOTM(NIY+2,i)*0.5*VYOTM(NIY,i),0.0);
end

%Upwind scheme for liquid x-axis momentum
VVYH=zeros(NIY+1,NIX+1);%H==half grid (cross grid)
for i=1:NIX+1
    for j=1:NIY+1
        VVYH(j,i)=VYOTM(j,i)*max(RFVXH(j,i),0.0)+VYOTM(j,i+1)*min(RFVXH(j,i),0.0);%including VVX3, VVX4 [kg/m/sec^2]
    end
end

%VX(1:NIY+2,1)=0.0 --> RFVXH(1:NIY+1,1)=0.0 --> VVYH(1:NIY+1,1)=0.0 --> Left impermeable boundary
%VX(1:NIY+2,NIX+1)=0.0 --> RFVXH(1:NIY+1,NIX+1)=0.0 --> VVYH(1:NIY+1,NIX+1)=0.0 --> Right impermeable boundary

%VY(1,1:NIX+2)=0.0 --> VVYH(1,1:NIX+1)=0.0 --> Top impermeable boundary
%VY(NIY+1,1:NIX+2)=0.0 --> VVYH(NIY+1,1:NIX+1)=0.0 --> Bottom impermeable boundary
%----------------------- 11. liquid y-flow momentum ---------------------------

%----------------------- 12. solid y-flow momentum ----------------------------

%----------------------- 12. solid y-flow momentum ----------------------------

%------------------------ 13. x-viscous momentum -----------------------------
%Originally marked as VDXI
VDSXI=zeros(NIY+2,NIX+2);%x velocity cell NORMAL STRESS [Pa.sec times 1/sec == Pa]
for i=2:NIX+1
    for j=1:NIY+2
        if(FSOTM(j,i)<=FSCR1)
            VDSXI(j,i)=MUOEI(j,i)*(VXOTM(j,i)-VXOTM(j,i-1))/dx(i-1);%including VDSX1, VDSX2
        else%Porus Region, crystals are interconnected
            %VDSXI(j,i)=MUOTM(j,i)*(VXOTM(j,i)*0.5*(FLOTM(j,i)+FLOTM(j,i+1))-VXOTM(j,i-1)*0.5*(FLOTM(j,i-1)+FLOTM(j,i)))/dx(i-1);%including VDX1, VDX2 of Xudaming
            VDSXI(j,i)=MUOTM(j,i)*FLOTM(j,i)*(VXOTM(j,i)-VXOTM(j,i-1))/dx(i-1);%Another expression of VDXI
        end
    end
end

for i=2:NIY+1
    if(FSOTM(i,1)<=FSCR1)%Crystal Free-Moving Region
        %IMPORTANT NOTE: since cavity is assumed rigid, any force NORMAL TO cavity wall from mixture is equal to the anti-force from cavity
        %wall, but with opposite direction (Newton's 3rd law), that is, VDSXI(1:NIY+2,1)==VDSXI(1:NIY+2,2) (for LEFT NORMAL STRESS) and
        %VDSXI(1:NIY+2,NIX+2)==VDSXI(1:NIY+2,NIX+1) (for RIGHT NORMAL STRESS). This is reduced to VXOTM(1:NIY+2,0)=-VXOTM(1:NIY+2,2).
        %Another explaination is: if VDSXI(1:NIY+2,1)<VDSXI(1:NIY+2,2), then some mass would be pushed out of the cavity through left wall.
        %So, to account for left impermeable boundary, we must set VXOTM(1:NIY+2,0)=-VXOTM(1:NIY+2,2).
        %We also set dx(0)=dx(1).
        VDSXI(i,1)=MUOEI(i,1)*(VXOTM(i,1)-(-VXOTM(i,2)))/dx(1);%left impermeable boundary

        %In the same manner, we set VXOTM(1:NIY+2,NIX+2)=-VXOTM(1:NIY+2,NIX) and set dx(NIX+1)=dx(NIX)
        VDSXI(i,NIX+2)=MUOEI(i,NIX+2)*(-VXOTM(i,NIX)-VXOTM(i,NIX+1))/dx(NIX);%right impermeable boundary
    else%Porus Region, crystals are interconnected
        %IMPORTANT NOTE: since cavity is assumed rigid, any force NORMAL TO cavity wall from liquid is equal to the anti-force from cavity
        %wall, but with opposite direction (Newton's 3rd law), that is, VDSXI(1:NIY+2,1)==VDSXI(1:NIY+2,2) (for LEFT NORMAL STRESS) and
        %VDSXI(1:NIY+2,NIX+2)==VDSXI(1:NIY+2,NIX+1) (for RIGHT NORMAL STRESS). This is reduced to VXOTM(1:NIY+2,0)=-VXOTM(1:NIY+2,2) and
        %FLOTM(1:NIY+2,0)=FLOTM(1:NIY+2,3). Another explaination is: if VDSXI(1:NIY+2,1)<VDSXI(1:NIY+2,2), then some mass would be pushed
        %out of the cavity through left wall. So, to account for left impermeable boundary, we must set VXOTM(1:NIY+2,0)=-VXOTM(1:NIY+2,2).
        %We also set dx(0)=dx(1).
        %VDSXI(i,1)=MUOTM(i,1)*(0.0-(-VXOTM(i,2)*0.5*(FLOTM(i,3)+FLOTM(i,2))))/dx(1);%left impermeable boundary
        %NOTE: FLOTM(i,3)+FLOTM(i,2) actually is FLOTM(i,0)+FLOTM(i,1)

        VDSXI(i,1)=MUOTM(i,1)*FLOTM(i,1)*(0.0-(-VXOTM(i,2)))/dx(1);%Another expression of VDXI

        %In the same manner, we set VXOTM(1:NIY+2,NIX+2)=-VXOTM(1:NIY+2,NIX) and FLOTM(1:NIY+2,NIX+3)=FLOTM(1:NIY+2,NIX). we also set dx(NIX)=dx(NIX+1).
        %VDSXI(i,NIX+2)=MUOTM(i,NIX+2)*(0.5*(FLOTM(i,NIX)+FLOTM(i,NIX+1))*(-VXOTM(i,NIX))-0.0)/dx(NIX);%right impermeable boundary
        %NOTE: FLOTM(i,NIX)+FLOTM(i,NIX+1) is actually FLOTM(i,NIX+2)+FLOTM(i,NIX+3)

        VDSXI(i,NIX+2)=MUOTM(i,NIX+2)*FLOTM(i,NIX+2)*((-VXOTM(i,NIX))-0.0)/dx(NIX);%Another expression of VDXI
    end
end

% %The following is a father set of the commented above
% MU=0.0;
% VXL=0.0;
% VXR=0.0;
% for i=2:NIX+1
%     for j=1:NIY+2
%         if(FSOTM(j,i)<=FSCR1)%Crystal Free-Moving Region
%             MU=MUOEI(j,i);
%         else%Porus Region, crystals are interconnected
%             MU=MUOTM(j,i);
%         end
%         if(0.5*(FSOTM(j,i-1)+FSOTM(j,i))<=FSCR1)
%             VXL=VXOTM(j,i-1);
%         else
%             VXL=0.5*(FLOTM(j,i-1)+FLOTM(j,i))*VXOTM(j,i-1);
%         end
%         if(0.5*(FSOTM(j,i)+FSOTM(j,i+1))<=FSCR1)
%             VXR=VXOTM(j,i);
%         else
%             VXR=0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i);
%         end
%         
%         VDSXI(j,i)=MU*(VXR-VXL)/dx(i-1);%including VDSX1, VDSX2
%     end
% end
% 
% %Left + Right impermeable boundary
% for i=2:NIY+1
%     %IMPORTANT NOTE: since cavity is assumed rigid, any force NORMAL TO cavity wall from mixture is equal to the anti-force from cavity
%     %wall, but with opposite direction (Newton's 3rd law), that is, VDSXI(1:NIY+2,1)==VDSXI(1:NIY+2,2) (for LEFT NORMAL STRESS) and
%     %VDSXI(1:NIY+2,NIX+2)==VDSXI(1:NIY+2,NIX+1) (for RIGHT NORMAL STRESS). This is reduced to VXOTM(1:NIY+2,0)=-VXOTM(1:NIY+2,2),
%     %FSOTM(1:NIY+2,0)=FSOTM(1:LNIY+2,3). Another explaination is: if VDSXI(1:NIY+2,1)<VDSXI(1:NIY+2,2), then some mass would be pushed out of the cavity through left wall.
%     %So, to account for left impermeable boundary, we must set VXOTM(1:NIY+2,0)=-VXOTM(1:NIY+2,2).
%     %We also set dx(0)=dx(1).
%     if(FSOTM(i,1)<=FSCR1)
%         MU=MUOEI(i,1);
%     else
%         MU=MUOTM(i,1);
%     end
%     if(0.5*(FSOTM(i,2)+FSOTM(i,3))<=FSCR1)
%         %FSOTM(i,3)+FSOTM(i,2) == FSOTM(i,0)+FSOTM(i,1)
%         VXL=-VXOTM(i,2);
%     else
%         VXL=-0.5*(FLOTM(i,2)+FLOTM(i,3))*VXOTM(i,2);
%         %FLOTM(i,3)+FLOTM(i,2) == FLOTM(i,0)+FLOTM(i,1)
%     end
%     
%     if(0.5*(FSOTM(i,1)+FSOTM(i,2))<=FSCR1)
%         VXR=VXOTM(i,1);
%     else
%         VXR=0.5*(FLOTM(i,1)+FLOTM(i,2))*VXOTM(i,1);
%     end
%     
%     VDSXI(i,1)=MU*(VXR-VXL)/dx(1);%VDSXI(i,1)=VDSXI(i,2) --> Left rigid/impermeable boundary
%     
%     %In the same manner, we set VXOTM(1:NIY+2,NIX+2)=-VXOTM(1:NIY+2,NIX) and FLOTM(1:NIY+2,NIX+3)=FLOTM(1:NIY+2,NIX). we also set dx(NIX)=dx(NIX+1).
%     if(FSOTM(i,NIX+2)<=FSCR1)
%         MU=MUOEI(i,NIX+2);
%     else
%         MU=MUOTM(i,NIX+2);
%     end
%     if(0.5*(FSOTM(i,NIX+1)+FSOTM(i,NIX+2))<=FSCR1)
%         VXL=VXOTM(i,NIX+1);
%     else
%         VXL=0.5*(FLOTM(i,NIX+1)+FLOTM(i,NIX+2))*VXOTM(i,NIX+1);
%     end
%     
%     if(0.5*(FSOTM(i,NIX)+FSOTM(i,NIX+1))<=FSCR1)
%         %FSOTM(i,NIX)+FSOTM(i,NIX+1) == FSOTM(i,NIX+2)+FSOTM(i,NIX+3)
%         VXR=-VXOTM(i,NIX);
%     else
%         VXR=-0.5*(FLOTM(i,NIX)+FLOTM(i,NIX+1))*VXOTM(i,NIX);
%         %FLOTM(i,NIX)+FLOTM(i,NIX+1) == FLOTM(i,NIX+2)+FLOTM(i,NIX+3)
%     end
%     
%     VDSXI(i,NIX+2)=MU*(VXR-VXL)/dx(NIX);%VDSXI(i,NIX+2)=VDSXI(i,NIX+1) --> Right rigid/impermeable boundary
%     
% end
% 
% %Left + Right impermeable boundary can be simply set as:
% %VDSXI(2:NIY+1,1)=VDSXI(2:NIY+1,2);
% %VDSXI(2:NIY+1,NIX+2)=VDSXI(2:NIY+1,NIX+1);

%Originally marked as VDXH
VDSXH=zeros(NIY+1,NIX+1);%x velocity cell SHEAR STRESS [Pa.sec times 1/sec == Pa]
for i=1:NIX+1
    for j=2:NIY

        if(0.25*(FSOTM(j,i)+FSOTM(j,i+1)+FSOTM(j+1,i)+FSOTM(j+1,i+1))<=FSCR1)%Free-Moving region
            VDSXH(j,i)=2.0*MUOEH(j,i)*(VXOTM(j+1,i)-VXOTM(j,i))/(dy(j-1)+dy(j));%[Pa.sec times 1/sec == Pa]

        else%Porus region, crystals are interconnected
            %VXOTM(j+1,i)~=0.0 and VXOTM(j,i)~=0.0
            %VXOTM(j+1,i)=0.0 and VXOTM(j,i)=0.0
            VDSXH(j,i)=2.0*MUOEH(j,i)*(0.5*(FLOTM(j+1,i)+FLOTM(j+1,i+1))*VXOTM(j+1,i)-0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i))/(dy(j-1)+dy(j));%[Pa.sec times 1/sec == Pa]

            %Up and Down row Logical for complete solid cells
            UDL1=logical((abs(VXOTM(j+1,i))<VErr)&&(abs(VXOTM(j,i))>=VErr));
            UDL2=logical(FLOTM(j+1,i+1)<FLErr||FLOTM(j+1,i)<FLErr);
            if(UDL1&&UDL2)
                %VXOTM(j+1,i)=0.0, VXOTM(j,i)~=0.0 and FLOTM(j+1,i+1)=0.0 or FLOTM(j+1,i)=0.0; set VXOTM(j+1,i)=-VXOTM(j,i)
                VDSXH(j,i)=2.0*MUOEH(j,i)*(0.5*(FLOTM(j+1,i)+FLOTM(j+1,i+1))*(-VXOTM(j,i))-0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i))/(dy(j-1)+dy(j));%[Pa.sec times 1/sec == Pa]
            end

            %Up and Down row Logical for complete solid cells
            UDL3=logical((abs(VXOTM(j+1,i))>=VErr)&&(abs(VXOTM(j,i))<VErr));
            UDL4=logical(FLOTM(j,i+1)<FLErr||FLOTM(j,i)<FLErr);
            if(UDL3&&UDL4)
                %VXOTM(j+1,i)~=0.0, VXOTM(j,i)=0.0 and FLOTM(j,i+1)=0.0 or FLOTM(j,i)=0.0; set VXOTM(j,i)=-VXOTM(j+1,i)
                VDSXH(j,i)=2.0*MUOEH(j,i)*(0.5*(FLOTM(j+1,i)+FLOTM(j+1,i+1))*(-VXOTM(j,i))-0.5*(FLOTM(j,i)+FLOTM(j,i+1))*(-VXOTM(j+1,i)))/(dy(j-1)+dy(j));%[Pa.sec times 1/sec == Pa]
            end
        end
    end
end

for i=1:NIX+1
    %VXOTM(1,i)=-VXOTM(2,i) top NO SLIP boundary
    VDSXH(1,i)=2.0*MUOEH(1,i)*(0.5*(FLOTM(2,i)+FLOTM(2,i+1))*VXOTM(2,i)-0.5*(FLOTM(1,i)+FLOTM(1,i+1))*VXOTM(1,i))/(2.0*dy(1));%VXOTM(1,i)=-VXOTM(2,i); [Pa.sec times 1/sec == Pa]
    
    %VXOTM(NIY+1,i)=VXOTM(NIY+2,i) bottom FREE boundary
    VDSXH(NIY+1,i)=2.0*MUOEH(NIY+1,i)*(0.5*(FLOTM(NIY+2,i)+FLOTM(NIY+2,i+1))*VXOTM(NIY+2,i)-0.5*(FLOTM(NIY+1,i)+FLOTM(NIY+1,i+1))*VXOTM(NIY+1,i))/(2.0*dy(NIY));%VXOTM(NIY+1,i)=-VXOTM(NIY+2,i);[Pa.sec times 1/sec == Pa]
end

% %The following is a father set of the commented above
% VXU=0.0;
% VXD=0.0;
% for i=1:NIX+1
%     for j=2:NIY
%         if(0.25*(FSOTM(j,i)+FSOTM(j,i+1)+FSOTM(j+1,i)+FSOTM(j+1,i+1))<=FSCR1)%Crystal Free-Moving Region
%             MU=MUOEH(j,i);
%         else%Porus Region, crystals are interconnected
%             MU=0.25*(MUOTM(j,i)+MUOTM(j,i+1)+MUOTM(j+1,i)+MUOTM(j+1,i+1));
%             %In theory, MU=MUOTM(j+0.5,i+0.5) calculated by VisCpRL(T(j+0.5,i+0.5),PA(j+0.5,i+0.5),MCL(j+0.5,i+0.5)). But it's too complex and may give similar results as presented here
%         end
%         if(0.5*(FSOTM(j,i)+FSOTM(j,i+1))<=FSCR1)
%             VXU=VXOTM(j,i);
%         else
%             VXU=0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i);
%         end
%         if(0.5*(FSOTM(j+1,i)+FSOTM(j+1,i+1))<=FSCR1)
%             VXD=VXOTM(j+1,i);
%         else
%             VXD=0.5*(FLOTM(j+1,i)+FLOTM(j+1,i+1))*VXOTM(j+1,i);
%         end
%         
%         VDSXH(j,i)=2.0*MU*(VXD-VXU)/(dy(j)+dy(j-1));%including VDSX1, VDSX2
%     end
% end
% 
% %Top + Bottom NO SLIP boundary: SOME SHEAR STRESS along x-axis wall
% for i=1:NIX+1
%     if(0.25*(FSOTM(1,i)+FSOTM(2,i)+FSOTM(1,i+1)+FSOTM(2,i+1))<=FSCR1)
%         MU=MUOEH(1,i);
%     else
%         MU=0.25*(MUOTM(1,i)+MUOTM(2,i)+MUOTM(1,i+1)+MUOTM(2,i+1));
%     end
%     if(0.5*(FSOTM(1,i)+FSOTM(1,i+1))<=FSCR1)
%         VXU=VXOTM(1,i);
%     else
%         VXU=0.5*(FLOTM(1,i)+FLOTM(1,i+1))*VXOTM(1,i);
%     end
%     if(0.5*(FSOTM(2,i)+FSOTM(2,i+1))<=FSCR1)
%         VXD=VXOTM(2,i);
%     else
%         VXD=0.5*(FLOTM(2,i)+FLOTM(2,i+1))*VXOTM(2,i);
%     end
%     
%     VDSXH(1,i)=2.0*MU*(VXD-VXU)/(dy(1)+dy(1));%VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> VDSXH(1,i)!=0.0 --> Top x-axis some shear stress
%     
%     if(0.25*(FSOTM(NIY+1,i)+FSOTM(NIY+2,i)+FSOTM(NIY+1,i+1)+FSOTM(NIY+2,i+1))<=FSCR1)
%         MU=MUOEH(NIY+1,i);
%     else
%         MU=0.25*(MUOTM(NIY+1,i)+MUOTM(NIY+2,i)+MUOTM(NIY+1,i+1)+MUOTM(NIY+2,i+1));
%     end
%     if(0.5*(FSOTM(NIY+1,i)+FSOTM(NIY+1,i+1))<=FSCR1)
%         VXU=VXOTM(NIY+1,i);
%     else
%         VXU=0.5*(FLOTM(NIY+1,i)+FLOTM(NIY+1,i+1))*VXOTM(NIY+1,i);
%     end
%     if(0.5*(FSOTM(NIY+2,i)+FSOTM(NIY+2,i+1))<=FSCR1)
%         VXD=VXOTM(NIY+2,i);
%     else
%         VXD=0.5*(FLOTM(NIY+2,i)+FLOTM(NIY+2,i+1))*VXOTM(NIY+2,i);
%     end
%     
%     VDSXH(NIY+1,i)=2.0*MU*(VXD-VXU)/(dy(NIY)+dy(NIY));%VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1)--> VDSXH(NIY+1,i)!=0.0 --> Bottom x-axis some shear stress
%     
% end
% 
% % %Top + Bottom FREE boundary: NO SHEAR STRESS along x-axis wall
% % for i=1:NIX+1
% %     if(0.25*(FSOTM(1,i)+FSOTM(2,i)+FSOTM(1,i+1)+FSOTM(2,i+1))<=FSCR1)
% %         MU=MUOEH(1,i);
% %     else
% %         MU=0.25*(MUOTM(1,i)+MUOTM(2,i)+MUOTM(1,i+1)+MUOTM(2,i+1));
% %     end
% %     if(0.5*(FSOTM(1,i)+FSOTM(1,i+1))<=FSCR1)
% %         VXU=VXOTM(1,i);
% %     else
% %         VXU=0.5*(FLOTM(1,i)+FLOTM(1,i+1))*VXOTM(1,i);
% %     end
% %     if(0.5*(FSOTM(2,i)+FSOTM(2,i+1))<=FSCR1)
% %         VXD=VXOTM(2,i);
% %     else
% %         VXD=0.5*(FLOTM(2,i)+FLOTM(2,i+1))*VXOTM(2,i);
% %     end
% %     
% %     VDSXH(1,i)=2.0*MU*(VXD-VXU)/(dy(1)+dy(1));%%VX(1,1:NIX+1)=VX(2,1:NIX+1) --> VDSXH(1,i)=0.0 --> Top x-axis NO shear stress
% %     
% %     if(0.25*(FSOTM(NIY+1,i)+FSOTM(NIY+2,i)+FSOTM(NIY+1,i+1)+FSOTM(NIY+2,i+1))<=FSCR1)
% %         MU=MUOEH(NIY+1,i);
% %     else
% %         MU=0.25*(MUOTM(NIY+1,i)+MUOTM(NIY+2,i)+MUOTM(NIY+1,i+1)+MUOTM(NIY+2,i+1));
% %     end
% %     if(0.5*(FSOTM(NIY+1,i)+FSOTM(NIY+1,i+1))<=FSCR1)
% %         VXU=VXOTM(NIY+1,i);
% %     else
% %         VXU=0.5*(FLOTM(NIY+1,i)+FLOTM(NIY+1,i+1))*VXOTM(NIY+1,i);
% %     end
% %     if(0.5*(FSOTM(NIY+2,i)+FSOTM(NIY+2,i+1))<=FSCR1)
% %         VXD=VXOTM(NIY+2,i);
% %     else
% %         VXD=0.5*(FLOTM(NIY+2,i)+FLOTM(NIY+2,i+1))*VXOTM(NIY+2,i);
% %     end
% %     
% %     VDSXH(NIY+1,i)=2.0*MU*(VXD-VXU)/(dy(NIY)+dy(NIY));%VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1)--> VDSXH(NIY+1,i)!=0.0 --> Bottom x-axis NO shear stress
% %     
% % end
%------------------------ 13. x-viscous momentum ----------------------------

%------------------------ 14. y-viscous momentum ----------------------------
%Originally marked as VDYI
VDSYI=zeros(NIY+2,NIX+2);%y velocity cell NORMAL STRESS [Pa.sec times 1/sec == Pa]
for i=1:NIX+2
    for j=2:NIY+1
        if(FSOTM(j,i)<=FSCR1)%Crystal Free-Moving Region
            VDSYI(j,i)=MUOEI(j,i)*(VYOTM(j,i)-VYOTM(j-1,i))/dy(j-1);
        else%Porus Region, crystals are interconnected
            VDSYI(j,i)=MUOTM(j,i)*(VYOTM(j,i)*0.5*(FLOTM(j+1,i)+FLOTM(j,i))-VYOTM(j-1,i)*0.5*(FLOTM(j,i)+FLOTM(j-1,i)))/dy(j-1);%including VDSY3, VDSY4
        end
    end
end

for i=2:NIX+1
    if(FSOTM(1,i)<=FSCR1)%Crystal Free-Moving Region
        %IMPORTANT NOTE: since cavity is assumed rigid, any force NORMAL TO cavity wall from mixture is equal to the anti-force from cavity
        %wall, but with opposite direction (Newton's 3rd law), that is, VDSYI(1,1:NIX+2)==VDSYI(2,1:NIX+2) (for TOP NORMAL STRSS) and
        %VDSYI(NIY+2,1:NIX+2)==VDSYI(NIY+1,1:NIX+2) (for BOTTOM NORMAL STRESS). This is reduced to VYOTM(0,1:NIX+2)=-VYOTM(2,1:NIX+2).
        %We also set dy(0)=dy(1).
        VDSYI(1,i)=MUOEI(1,i)*(0.0-(-VYOTM(2,i)))/dy(1);%top impermeable boundary

        %In the same manner, we set VYOTM(NIY+2,1:NIX+2)=-VYOTM(NIY,1:NIX+2) and set dy(NIY+1)=dy(NIY).
        VDSYI(NIY+2,i)=MUOEI(NIY+2,i)*((-VYOTM(NIY,i))-0.0)/dy(NIY);%bottom impermeable boundary

    else%Porus region, crystals are interconnected
        %IMPORTANT NOTE: since cavity is assumed rigid, any force NORMAL TO cavity wall from mixture is equal to the anti-force from cavity
        %wall, but with opposite direction (Newton's 3rd law), that is, VDSYI(1,1:NIX+2)==VDSYI(2,1:NIX+2) (for TOP NORMAL STRESS) and
        %VDSYI(NIY+2,1:NIX+2)==VDSYI(NIY+1,1:NIX+2) (for BOTTOM NORMAL STRESS). This is reduced to VYOTM(0,1:NIX+2)=-VYOTM(2,1:NIX+2) and
        %FLOTM(0,1:NIX+2)=FLOTM(3,1:NIX+2). We also set dy(0)=dy(1).
        VDSYI(1,i)=MUOTM(1,i)*(0.0-0.5*(FLOTM(2,i)+FLOTM(3,i))*(-VYOTM(2,i)))/dy(1);%top impermeable boundary
        %NOTE: FLOTM(2,i)+FLOTM(3,i) is actually FLOTM(1,i)+FLOTM(0,i)

        %In the same manner, we set VYOTM(NIY+2,1:NIX+2)=-VYOTM(NIY,1:NIX+2) and FLOTM(NIY+3,1:NIX+2)=FLOTM(NIY,1:NIX+2), we also set dy(NIY)=dy(NIY+1).
        VDSYI(NIY+2,i)=MUOTM(NIY+2,i)*(0.5*(FLOTM(NIY+1,i)+FLOTM(NIY,i))*(-VYOTM(NIY,i))-0.0)/dy(NIY);%bottom impermeable boundary
        %NOTE: FLOTM(NIY+1,i)+FLOTM(NIY,i) is actually FLOTM(NIY+3,i)+FLOTM(NIY+2,i)
    end
end

% VYU=0.0;
% VYD=0.0;
% for i=1:NIX+2
%     for j=2:NIY+1
%         if(FSOTM(j,i)<=FSCR1)
%             MU=MUOEI(j,i);
%         else
%             MU=MUOTM(j,i);
%         end
%         if(0.5*(FSOTM(j-1,i)+FSOTM(j,i))<=FSCR1)
%             VYU=VYOTM(j-1,i);
%         else
%             VYU=0.5*(FLOTM(j-1,i)+FLOTM(j,i))*VYOTM(j-1,i);
%         end
%         if(0.5*(FSOTM(j,i)+FSOTM(j+1,i))<=FSCR1)
%             VYD=VYOTM(j,i);
%         else
%             VYD=0.5*(FLOTM(j,i)+FLOTM(j+1,i))*VYOTM(j,i);
%         end
%         
%         VDSYI(j,i)=MU*(VYD-VYU)/dy(j-1);
%     end
% end
% 
% %Top + Bottom impermeable boundary
% for i=2:NIX+1
%     if(FSOTM(1,i)<=FSCR1)
%         MU=MUOEI(1,i);
%     else
%         MU=MUOTM(1,i);
%     end
%     if(0.5*(FSOTM(2,i)+FSOTM(3,i))<=FSCR1)
%         %FSOTM(2,i)+FSOTM(3,i) == FSOTM(1,i)+FSOTM(0,i)
%         VYU=-VYOTM(2,i);
%     else
%         VYU=-0.5*(FLOTM(2,i)+FLOTM(3,i))*VYOTM(2,i);
%         %FLOTM(2,i)+FLOTM(3,i) == FLOTM(1,i)+FLOTM(0,i)
%     end
%     if(0.5*(FSOTM(1,i)+FSOTM(2,i))<=FSCR1)
%         VYD=VYOTM(1,i);
%     else
%         VYD=0.5*(FLOTM(1,i)+FLOTM(2,i))*VYOTM(1,i);
%     end
%     
%     VDSYI(1,i)=MU*(VYD-VYU)/dy(1);%VDSYI(1,i)=VDSYI(2,i) --> Top rigid/impermeable boundary 
%     
%     if(FSOTM(NIY+2,i)<=FSCR1)
%         MU=MUOEI(NIY+2,i);
%     else
%         MU=MUOTM(NIY+2,i);
%     end
%     if(0.5*(FSOTM(NIY+1,i)+FSOTM(NIY+2,i))<=FSCR1)
%         VYU=VYOTM(NIY+1,i);
%     else
%         VYU=0.5*(FLOTM(NIY+1,i)+FLOTM(NIY+2,i))*VYOTM(NIY+1,i);
%     end
%     if(0.5*(FSOTM(NIY,i)+FSOTM(NIY+1,i))<=FSCR1)
%         %FSOTM(NIY+1,i)+FSOTM(NIY,i) == FSOTM(NIY+3,i)+FSOTM(NIY+2,i)
%         VYD=-VYOTM(NIY,i);
%     else
%         VYD=-0.5*(FLOTM(NIY,i)+FLOTM(NIY+1,i))*VYOTM(NIY,i);
%         %FLOTM(NIY+1,i)+FLOTM(NIY,i) == FLOTM(NIY+3,i)+FLOTM(NIY+2,i)
%     end
%     
%     VDSYI(NIY+2,i)=MU*(VYD-VYU)/dy(NIY);%VDSYI(NIY+2,i)=VDSYI(NIY+1,i) --> Bottom rigid/impermeable boundary
%     
% end

VDSYH=zeros(NIY+1,NIX+1);%y velocity cell SHEAR STRESS [Pa.sec times 1/sec == Pa]
for i=2:NIX
    %VYOTM(1,i)==0.0 --> top impermeable boundary --> VDSYH(1,i)=0.0
    %VYOTM(NIY+1,i)==0.0 --> bottom impermeable boundary --> VDSYH(NIY+1,i)=0.0
    for j=1:NIY+1
        %VYOTM(j,i)~=0.0 and VYOTM(j,i+1)~=0.0
        %VYOTM(j,i)=0.0 and VYOTM(j,i+1)=0.0
        VDSYH(j,i)=2.0*0.25*(MUOTM(j,i)+MUOTM(j,i+1)+MUOTM(j+1,i)+MUOTM(j+1,i+1))*(0.5*(FLOTM(j+1,i+1)+FLOTM(j,i+1))*VYOTM(j,i+1)-0.5*(FLOTM(j+1,i)+FLOTM(j,i))*VYOTM(j,i))/(dx(i-1)+dx(i));%[Pa.sec times 1/sec == Pa]

        %Left and Right column Logical
        LRL1=logical(abs(VYOTM(j,i))<VErr&&abs(VYOTM(j,i+1))>=VErr);
        LRL2=logical(FLOTM(j,i)<FLErr||FLOTM(j+1,i)<FLErr);
        if(LRL1&&LRL2)
            %VYOTM(j,i)=0.0 and VYOTM(j,i+1)~=0.0 and FLOTM(j,i)=0.0 or FLOTM(j+1,i)=0.0; set VYOTM(j,i)=-VYOTM(j,i+1)
            VDSYH(j,i)=2.0*0.25*(MUOTM(j,i)+MUOTM(j,i+1)+MUOTM(j+1,i)+MUOTM(j+1,i+1))*(0.5*(FLOTM(j+1,i+1)+FLOTM(j,i+1))*VYOTM(j,i+1)-0.5*(FLOTM(j+1,i)+FLOTM(j,i))*(-VYOTM(j,i+1)))/(dx(i-1)+dx(i));%[Pa.sec times 1/sec == Pa]
        end

        %Left and Right column Logical
        LRL3=logical(abs(VYOTM(j,i))>=VErr&&abs(VYOTM(j,i+1))<VErr);
        LRL4=logical(FLOTM(j,i+1)<FLErr||FLOTM(j+1,i+1)<FLErr);
        if(LRL3&&LRL4)
            %VYOTM(j,i)~=0.0 and VYOTM(j,i+1)=0.0; set VYOTM(j,i+1)=-VYOTM(j,i)
            VDSYH(j,i)=2.0*0.25*(MUOTM(j,i)+MUOTM(j,i+1)+MUOTM(j+1,i)+MUOTM(j+1,i+1))*(0.5*(FLOTM(j+1,i+1)+FLOTM(j,i+1))*(-VYOTM(j,i))-0.5*(FLOTM(j+1,i)+FLOTM(j,i))*(-VYOTM(j,i+1)))/(dx(i-1)+dx(i));%[Pa.sec times 1/sec == Pa]
        end

    end
end

for i=1:NIY+1
    %VYOTM(i,1)==VYOTM(i,2) --> left FREE boundary
    VDSYH(i,1)=2.0*0.25*(MUOTM(i,1)+MUOTM(i,2)+MUOTM(i+1,1)+MUOTM(i+1,2))*(0.5*(FLOTM(i+1,2)+FLOTM(i,2))*VYOTM(i,2)-0.5*(FLOTM(i+1,1)+FLOTM(i,1))*VYOTM(i,1))/(2.0*dx(1));%VYOTM(i,1)=VYOTM(i,2) for slip boundary; [Pa.sec times 1/sec == Pa]
    
    %VYOTM(i,NIX+2)=VYOTM(i,NIX+2) --> right FREE boundary
    VDSYH(i,NIX+1)=2.0*0.25*(MUOTM(i,NIX+1)+MUOTM(i,NIX+2)+MUOTM(i+1,NIX+1)+MUOTM(i+1,NIX+2))*(0.5*(FLOTM(i+1,NIX+2)+FLOTM(i,NIX+2))*VYOTM(i,NIX+2)-0.5*(FLOTM(i+1,NIX+1)+FLOTM(i,NIX+1))*VYOTM(i,NIX+1))/(2.0*dx(NIX));%VYOTM(i,NIX+2)=VYOTM(i,NIX+1); [Pa.sec times 1/sec == Pa]
end

% VYL=0.0;
% VYR=0.0;
% for i=2:NIX
%     for j=1:NIY+1
%         if(0.25*(FSOTM(j,i)+FSOTM(j+1,i)+FSOTM(j,i+1)+FSOTM(j+1,i+1))<=FSCR1)
%             MU=MUOEH(j,i);
%         else
%             MU=0.25*(MUOTM(j,i)+MUOTM(j+1,i)+MUOTM(j,i+1)+MUOTM(j+1,i+1));
%         end
%         if(0.5*(FSOTM(j,i)+FSOTM(j+1,i))<=FSCR1)
%             VYL=VYOTM(j,i);
%         else
%             VYL=0.5*(FLOTM(j,i)+FLOTM(j+1,i))*VYOTM(j,i);
%         end
%         if(0.5*(FSOTM(j,i+1)+FSOTM(j+1,i+1))<=FSCR1)
%             VYR=VYOTM(j,i+1);
%         else
%             VYR=0.5*(FLOTM(j,i+1)+FLOTM(j+1,i+1))*VYOTM(j,i+1);
%         end
%         VDSYH(j,i)=2.0*MU*(VYR-VYL)/(dx(i-1)+dx(i));
%     end
% end
% 
% %Left + Right FREE boundary: NO SHEAR STRESS along y-axis wall, that is, VY(1:NIY+1,1)=VY(1:NIY+1,2) or VY(1:NIY+1,NIX+2)=VY(1:NIY+1,NIX+1)
% for i=1:NIY+1
%     if(0.25*(FSOTM(i,1)+FSOTM(i+1,1)+FSOTM(i,2)+FSOTM(i+1,2))<=FSCR1)
%         MU=MUOEH(i,1);
%     else
%         MU=0.25*(MUOTM(i,1)+MUOTM(i+1,1)+MUOTM(i,2)+MUOTM(i+1,2));
%     end
%     if(0.5*(FSOTM(i,1)+FSOTM(i+1,1))<=FSCR1)
%         VYL=VYOTM(i,1);
%     else
%         VYL=0.5*(FLOTM(i,1)+FLOTM(i+1,1))*VYOTM(i,1);
%     end
%     if(0.5*(FSOTM(i,2)+FSOTM(i+1,2))<=FSCR1)
%         VYR=VYOTM(i,2);
%     else
%         VYR=0.5*(FLOTM(i,2)+FLOTM(i+1,2))*VYOTM(i,2);
%     end
%     VDSYH(i,1)=2.0*MU*(VYR-VYL)/(dx(1)+dx(1));%VDSYH(i,1)=0.0 --> NO SHEAR STRESS along y-axis wall --> Left FREE boundary
%     
%     if(0.25*(FSOTM(i,NIX+1)+FSOTM(i+1,NIX+1)+FSOTM(i,NIX+2)+FSOTM(i+1,NIX+2))<=FSCR1)
%         MU=MUOEH(i,NIX+1);
%     else
%         MU=0.25*(MUOTM(i,NIX+1)+MUOTM(i+1,NIX+1)+MUOTM(i,NIX+2)+MUOTM(i+1,NIX+2));
%     end
%     if(0.5*(FSOTM(i,NIX+1)+FSOTM(i+1,NIX+1))<=FSCR1)
%         VYL=VYOTM(i,NIX+1);
%     else
%         VYL=0.5*(FLOTM(i,NIX+1)+FLOTM(i+1,NIX+1))*VYOTM(i,NIX+1);
%     end
%     if(0.5*(FSOTM(i,NIX+2)+FSOTM(i+1,NIX+2))<=FSCR1)
%         VYR=VYOTM(i,NIX+2);
%     else
%         VYR=0.5*(FLOTM(i,NIX+2)+FLOTM(i+1,NIX+2))*VYOTM(i,NIX+2);
%     end
%     VDSYH(i,NIX+1)=2.0*MU*(VYR-VYL)/(dx(NIX)+dx(NIX));%VDSYH(i,NIX+1)=0.0 --> NO SHEAR STRESS along y-axis wall --> Right FREE boundary
%     
% end
% 
% % %Left + Right NO SLIP boundary: SOME SHEAR STRESS along y-axis wall, that is, VY(1:NIY+1,1)=-VY(1:NIY+1,2) or VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1)
% % for i=1:NIY+1
% %     if(0.25*(FSOTM(i,1)+FSOTM(i+1,1)+FSOTM(i,2)+FSOTM(i+1,2))<=FSCR1)
% %         MU=MUOEH(j,1);
% %     else
% %         MU=0.25*(MUOTM(i,1)+MUOTM(i+1,1)+MUOTM(i,2)+MUOTM(i+1,2));
% %     end
% %     if(0.5*(FSOTM(i,1)+FSOTM(i+1,1))<=FSCR1)
% %         VYL=VYOTM(i,1);
% %     else
% %         VYL=0.5*(FLOTM(i,1)+FLOTM(i+1,1))*VYOTM(i,1);
% %     end
% %     if(0.5*(FSOTM(i,2)+FSOTM(i+1,2))<=FSCR1)
% %         VYR=VYOTM(i,2);
% %     else
% %         VYR=0.5*(FLOTM(i,2)+FLOTM(i+1,2))*VYOTM(i,2);
% %     end
% %     VDSYH(i,1)=2.0*MU*(VYR-VYL)/(dx(1)+dx(1));%VDSYH(i,1)!=0.0 --> Some SHEAR STRESS along y-axis wall --> Left NO SLIP boundary
% %     
% %     if(0.25*(FSOTM(i,NIX+1)+FSOTM(i+1,NIX+1)+FSOTM(i,NIX+2)+FSOTM(i+1,NIX+2))<=FSCR1)
% %         MU=MUOEH(j,NIX+1);
% %     else
% %         MU=0.25*(MUOTM(i,NIX+1)+MUOTM(i+1,NIX+1)+MUOTM(i,NIX+2)+MUOTM(i+1,NIX+2));
% %     end
% %     if(0.5*(FSOTM(i,NIX+1)+FSOTM(i+1,NIX+1))<=FSCR1)
% %         VYL=VYOTM(i,NIX+1);
% %     else
% %         VYL=0.5*(FLOTM(i,NIX+1)+FLOTM(i+1,NIX+1))*VYOTM(i,NIX+1);
% %     end
% %     if(0.5*(FSOTM(i,NIX+2)+FSOTM(i+1,NIX+2))<=FSCR1)
% %         VYR=VYOTM(i,NIX+2);
% %     else
% %         VYR=0.5*(FLOTM(i,NIX+2)+FLOTM(i+1,NIX+2))*VYOTM(i,NIX+2);
% %     end
% %     VDSYH(i,NIX+1)=2.0*MU*(VYR-VYL)/(dx(NIX)+dx(NIX));%VDSYH(i,NIX+1)!=0.0 --> Some SHEAR STRESS along y-axis wall --> Right NO SLIP boundary
% %     
% % end
%------------------------ 14. y-viscous momentum ----------------------------

%----------------------- 15. x-convection liquid ----------------------------
VXCONV=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        VXCONV(j,i)=-2.0*dtb*(VVXI(j+1,i+1)-VVXI(j+1,i))/(dx(i)+dx(i-1))-dtb*(VVXH(j+1,i)-VVXH(j,i))/dy(j);%[kg/m^2/sec]
        
        %For solid cells, their faces are not permeable, thus there are not
        %convective flux through these faces.
        %NOTE: We can not set VXCONV(j,i)=0 by setting
        %VVXI(j+1,i+1)==VVXI(j+1,i) because VVXI needs RFVXI which may have
        %two different values at the same time for mirror projection to
        %ensure the faces have 0.0 convective flux.
        %Left and Right column Logical
        LRL1=logical(abs(VXOTM(j+1,i))<VErr);
        LRL2=logical(FLOTM(j+1,i)<FLErr||FLOTM(j+1,i+1)<FLErr);
        if(LRL1&&LRL2)
            VXCONV(j,i)=0.0;
        end
        
    end
end

%Recheck boundary condition explicitly
for i=1:NIY
    %VVXI(1:NIY+2,1)==VVXI(1:NIY+2,2), VVXH(1:NIY+1,1)=0.0 --> VXCONV(i,1)=0.0 --> Left impermeable boundary 
    VXCONV(i,1)=-2.0*dtb*(VVXI(i+1,2)-VVXI(i+1,1))/(2.0*dx(1))-dtb*(VVXH(i+1,1)-VVXH(i,1))/dy(i);%[kg/m^2/sec]
    
    %VVXI(1:NIY+2,NIX+2)==VVXI(1:NIY+2,NIX+1), VVXH(1:NIY+1,NIX+1)=0.0 --> VXCONV(i,NIX+1)=0.0 --> Right impermeable boundary 
    VXCONV(i,NIX+1)=-2.0*dtb*(VVXI(i+1,NIX+2)-VVXI(i+1,NIX+1))/(2.0*dx(NIX))-dtb*(VVXH(i+1,NIX+1)-VVXH(i,NIX+1))/dy(i);%[kg/m^2/sec]
end
%----------------------- 15. x-convection liquid ----------------------------

%----------------------- 16. x-convection solid -----------------------------

%----------------------- 16. x-convection solid -----------------------------

%----------------------- 17. x-diffusion mixture ----------------------------
%Originally marked as VXDIF
VSXDIF=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        VSXDIF(j,i)=2.0*dtb*(VDSXI(j+1,i+1)-VDSXI(j+1,i))/(dx(i)+dx(i-1))+dtb*(VDSXH(j+1,i)-VDSXH(j,i))/dy(j);%[Pa.sec/m]
        
        %For solid cells, their faces are not permeable, thus there are not
        %diffusive flux through these faces.
        %NOTE: We can not set VSXDIF(j,i)=0 by setting
        %VDSXI(j+1,i+1)==VDSXI(j+1,i) and VDSXH(j+1,i)==VDSXH(j,i) because VDSXI
        %and VDSXH are products of primitive varibles, thus VDSXI and VDSXH
        %are in fact right, so we have to set VSXDIF manually.
        %Left and Right column Logical
        LRL1=logical(abs(VXOTM(j+1,i))<VErr);
        LRL2=logical(FLOTM(j+1,i)<FLErr||FLOTM(j+1,i+1)<FLErr);
        if(LRL1&&LRL2)
            VSXDIF(j,i)=0.0;
        end
        
    end
end

%Recheck boundary condition explicitly
for i=1:NIY
    %VDSXI(1:NIY+2,1)=VDSXI(1:NIY+2,2) --> NORMAL STRESS NORMAL TO left wall --> Left rigid/impermeable boundary
    %VX(1:NIY+2,1)=0.0 --> VDSXH(1:NIY+1,1)=0.0 --> Left impermeable boundary --> VSXDIF(1:NIY,1)=0.0
    VSXDIF(i,1)=2.0*dtb*(VDSXI(i+1,2)-VDSXI(i+1,1))/(2.0*dx(1))+dtb*(VDSXH(i+1,1)-VDSXH(i,1))/dy(i);%[Pa.sec/m]
    
    %VDSXI(1:NIY+2,NIX+2)=VDSXI(1:NIY+2,NIX+1) --> NORMAL STRESS NORMAL TO right wall --> Right rigid/impermeable boundary
    %VX(1:NIY+2,NIX+1)=0.0 --> VDSXH(1:NIY+1,NIX+1)=0.0 --> Right impermeable boundary --> VSXDIF(1:NIY,NIX+1)=0.0
    VSXDIF(i,NIX+1)=2.0*dtb*(VDSXI(i+1,NIX+2)-VDSXI(i+1,NIX+1))/(2.0*dx(NIX))+dtb*(VDSXH(i+1,NIX+1)-VDSXH(i,NIX+1))/dy(i);%[Pa.sec/m]
    
end
%----------------------- 17. x-diffusion mixture ----------------------------

%----------------------- 18. y-convection liquid ----------------------------
VYCONV=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        VYCONV(j,i)=-dtb*(VVYH(j,i+1)-VVYH(j,i))/dx(i)-2.0*dtb*(VVYI(j+1,i+1)-VVYI(j,i+1))/(dy(j)+dy(j-1));%[Pa.sec/m]
        
        %For solid cells, their faces are not permeable, thus there are not
        %convective flux through these faces.
        %NOTE: We can not set VYCONV(j,i)=0 by setting
        %VVYI(j+1,i+1)==VVYI(j,i+1) because VVYI needs RFVYI which may have
        %two different values at the same time for mirror projection to
        %ensure the faces have 0.0 convective flux.
        %Up and Down column Logical
        UDL1=logical(abs(VYOTM(j,i+1))<VErr);
        UDL2=logical(FLOTM(j,i+1)<FLErr||FLOTM(j+1,i+1)<FLErr);
        if(UDL1&&UDL2)
            VYCONV(j,i)=0.0;
        end
        
    end
end

%Recheck boundary condition explicitly
for i=1:NIX
    %VY(1,1:NIX+2)=0.0 --> VVYH(1,1:NIX+1)=0.0 --> Top impermeable boundary 
    %VVYI(1,1:NIX+2)=VVYI(2,1:NIX+2) --> VYCONV(1,1:NIX)=0.0
    VYCONV(1,i)=-dtb*(VVYH(1,i+1)-VVYH(1,i))/dx(i)-2.0*dtb*(VVYI(2,i+1)-VVYI(1,i+1))/(2.0*dy(1));%[Pa.sec/m]
    
    %VY(NIY+1,1:NIX+2)=0.0 --> VVYH(NIY+1,1:NIX+1)=0.0 --> Bottom impermeable boundary
    %VVYI(NIY+2,1:NIX+2)=VVYI(NIY+1,1:NIX+2) --> VYCONV(NIY+1,1:NIX)=0.0
    VYCONV(NIY+1,i)=-dtb*(VVYH(NIY+1,i+1)-VVYH(NIY+1,i))/dx(i)-2.0*dtb*(VVYI(NIY+2,i+1)-VVYI(NIY+1,i+1))/(2.0*dy(NIY));%[Pa.sec/m]
end
%----------------------- 18. y-convection liquid ----------------------------

%----------------------- 19. y-convection solid -----------------------------

%----------------------- 19. y-convection solid -----------------------------

%----------------------- 20. y-diffusion mixture ----------------------------
%Originally marked as VYDIF
VSYDIF=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        VSYDIF(j,i)=dtb*(VDSYH(j,i+1)-VDSYH(j,i))/dx(i)+2.0*dtb*(VDSYI(j+1,i+1)-VDSYI(j,i+1))/(dy(j)+dy(j-1));%[Pa.sec/m]
        
        %For solid cells, their faces are not permeable, thus there are not
        %diffusive flux through these faces.
        %NOTE: We can not set VSYDIF(j,i)=0 by setting
        %VDSYI(j+1,i+1)==VDSYI(j,i+1) and VDSYH(j,i+1)==VDSYH(j,i) because VDSYI
        %and VDSYH are products of primitive varibles, thus VDSYI and VDSYH
        %are in fact right, so we have to set VSYDIF manually
        %Up and Down column Logical
        UDL1=logical(abs(VYOTM(j,i+1))<VErr);
        UDL2=logical(FLOTM(j,i+1)<FLErr||FLOTM(j+1,i+1)<FLErr);
        if(UDL1&&UDL2)
            VSYDIF(j,i)=0.0;
        end
        
    end
end

%Recheck boundary condition explicitly
for i=1:NIX
    %VY(1,1:NIX+2)=0.0 --> VDSYH(1,1:NIX+1)=0.0 --> Top impermeable bounadry
    %VDSYI(2,1:NIX+2)=VDSYI(1,1;NIX+2) --> NORMAL STRESS NORMAL TO top wall --> VSYDIF(1,1:NIX)=0.0
    VSYDIF(1,i)=dtb*(VDSYH(1,i+1)-VDSYH(1,i))/dx(i)+2.0*dtb*(VDSYI(2,i+1)-VDSYI(1,i+1))/(2.0*dy(1));%[Pa.sec/m]
    
    %VY(NIY+1,1:NIX+2)=0.0 --> VDSYH(NIY+1,1:NIX+1)=0.0 --> Bottom impermeable boundary
    %VDSYI(NIY+2,1:NIX+2)=VDSYI(NIY+1,1:NIX+2) --> NORMAL STRESS NORMAL TO bottom wall --> VSYDIF(NIY+1,1:NIX)=0.0
    VSYDIF(NIY+1,i)=dtb*(VDSYH(NIY+1,i+1)-VDSYH(NIY+1,i))/dx(i)+2.0*dtb*(VDSYI(NIY+2,i+1)-VDSYI(NIY+1,i+1))/(2.0*dy(NIY));%[Pa.sec/m]
end
%----------------------- 20. y-diffusion mixture ----------------------------

%----------------------- 21. x- pressure gradient ---------------------------
FPDXL=0.0;
FPDXR=0.0;
FPDX=zeros(NIY+2,NIX+1);%D==Difference
%NOTE: FPDX will be used in RFVXTM, while RFVXTM uses old pressure to guess (see PPT P.64)
for i=2:NIX
    for j=1:NIY+2
        if(FSOTM(j,i)<=FSCR1)
            FPDXL=PDTM(j,i);
        else
            FPDXL=FLOTM(j,i)*PDTM(j,i);
            %If FLOTM(j,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        if(FSOTM(j,i+1)<=FSCR1)
            FPDXR=PDTM(j,i+1);
        else
            FPDXR=FLOTM(j,i+1)*PDTM(j,i+1);
            %If FLOTM(j,i+1)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        
        FPDX(j,i)=-2.0*dtb*(FPDXR-FPDXL)/(dx(i)+dx(i-1));%[Pa.sec/m == kg/m^2/sec]
        
        %solid cell at (j,i), set FPDX=0.0 mandatorily
        if(FLOTM(j,i)<FLErr)
            FPDX(j,i-1)=0.0;%cell left face
            FPDX(j,i)=0.0;%cell right face
        end
        
    end
end

%Left & Right boundary
for i=1:NIY+2
    if(FSOTM(i,1)<=FSCR1)
        FPDXL=PDTM(i,1);
    else
        FPDXL=FLOTM(i,1)*PDTM(i,1);
        %If FLOTM(j,1)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSOTM(i,2)<=FSCR1)
        FPDXR=PDTM(i,2);
    else
        FPDXR=FLOTM(i,2)*PDTM(i,2);
        %If FLOTM(j,2)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    FPDX(i,1)=-2.0*dtb*(FPDXR-FPDXL)/(2.0*dx(1));%left boundary [Pa.sec/m == kg/m^2/sec]
    %NOTE: FLOTM(1:NIY+2,1)=FLOTM(1:NIY+2,2),
    %PDTM(1:NIY+2,1)=PDTM(1:NIY+2,2) --> FPDX(1:NIY+2,1)=0.0
    
    %solid at cell (i,2), set FPDX=0.0 mandatorily
    if(FLOTM(i,2)<=FLErr)
        FPDX(i,1)=0.0;%cell left face
        FPDX(i,2)=0.0;%cell right face
    end
    
    if(FSOTM(i,NIX+1)<=FSCR1)
        FPDXL=PDTM(i,NIX+1);
    else
        FPDXL=FLOTM(i,NIX+1)*PDTM(i,NIX+1);
        %If FLOTM(j,NIX+1)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSOTM(i,NIX+2)<=FSCR1)
        FPDXR=PDTM(i,NIX+2);
    else
        FPDXR=FLOTM(i,NIX+2)*PDTM(i,NIX+2);
        %If FLOTM(j,NIX+2)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    FPDX(i,NIX+1)=-2.0*dtb*(FPDXR-FPDXL)/(2.0*dx(NIX));%right boundary [Pa.sec/m == kg/m^2/sec]
    %NOTE: FLOTM(1:NIY+2,NIX+1)=FLOTM(1:NIY+2,NIX+2),
    %PDTM(1:NIY+2,NIX+1)=PDTM(1:NIY+2,NIX+2) --> FPDX(1:NIY+2,NIX+1)=0.0

    %solid at cell (i,NIX+1), set FPDX=0.0 mandatorily
    if(FLOTM(i,NIX+1)<=FLErr)
        FPDX(i,NIX)=0.0;%cell left face
        FPDX(i,NIX+1)=0.0;%cell right face
    end
    
end

FP0XL=0.0;
FP0XR=0.0;
PSFX=zeros(NIY+2,NIX+1);%Initial static pressure * FLN, to extract static pressure from actual pressure
%NOTE: FL is FLN (New step FL)!
for i=2:NIX
    for j=1:NIY+2
        if(FSNTM(j,i)<=FSCR1)
            FP0XL=P0(j,i);
        else
            FP0XL=FLNTM(j,i)*P0(j,i);
            %If FLNTM(j,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        if(FSNTM(j,i+1)<=FSCR1)
            FP0XR=P0(j,i+1);
        else
            FP0XR=FLNTM(j,i+1)*P0(j,i+1);
            %If FLNTM(j,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        
        PSFX(j,i)=-2.0*dtb*(FP0XR-FP0XL)/(dx(i)+dx(i-1));%[Pa.sec/m == kg/m^2/sec]
        
        %If VXOTM(j,i)=0.0, effective static pressure does not work, which
        %means that not any flux can go through cell face
        
        if(FLOTM(j,i)<FLErr)
            PSFX(j,i-1)=0.0;
            PSFX(j,i)=0.0;
        end
        
    end
end

%Left & Right boundary
for i=1:NIY+2
    if(FSNTM(i,1)<=FSCR1)
        FP0XL=P0(i,1);
    else
        FP0XL=FLNTM(i,1)*P0(i,1);
        %If FLNTM(j,1)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSNTM(i,2)<=FSCR1)
        FP0XR=P0(i,2);
    else
        FP0XR=FLNTM(i,2)*P0(i,2);
        %If FLNTM(j,2)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    PSFX(i,1)=-2.0*dtb*(FP0XR-FP0XL)/(2.0*dx(1));%left boundary [Pa.sec/m == kg/m^2/sec]
    %P0(1:NIY+2,1)=P0(1:NIY+2,2), FLNTM(1:NIY+2,1)=FLNTM(1:NIY+2,2) --> PSFX(1:NIY+2,1)=0.0
    
    %solid at cell (i,2), set PSFX=0.0 mandatorily
    if(FLOTM(i,2)<FLErr)
        PSFX(i,1)=0.0;
        PSFX(i,2)=0.0;
    end
    
    if(FSNTM(i,NIX+1)<=FSCR1)
        FP0XL=P0(i,NIX+1);
    else
        FP0XL=FLNTM(i,NIX+1)*P0(i,NIX+1);
        %If FLNTM(j,NIX+1)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSNTM(i,NIX+2)<=FSCR1)
        FP0XR=P0(i,NIX+2);
    else
        FP0XR=FLNTM(i,NIX+2)*P0(i,NIX+2);
        %If FLNTM(j,NIX+2)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    PSFX(i,NIX+1)=-2.0*dtb*(FP0XR-FP0XL)/(2.0*dx(NIX));%right boundary [Pa.sec/m == kg/m^2/sec]
    %P0(1:NIY+2,NIX+1)=P0(1:NIY+2,NIX+2), FLNTM(1:NIY+2,NIX+1)=FLNTM(1:NIY+2,NIX+2) --> PSFX(1:NIY+2,NIX+1)=0.0

    %solid at cell (i,NIX+1), set PSFX=0.0 mandatorily
    if(FLOTM(i,NIX+1)<FLErr)
        PSFX(i,NIX)=0.0;
        PSFX(i,NIX+1)=0.0;
    end
    
end
%----------------------- 21. x- pressure gradient ---------------------------

%----------------------- 22. y- pressure gradient ---------------------------
FPDYU=0.0;
FPDYD=0.0;
FPDY=zeros(NIY+1,NIX+2);%D==Difference
for i=1:NIX+2
    for j=2:NIY
        if(FSOTM(j,i)<=FSCR1)
            FPDYU=PDTM(j,i);
        else
            FPDYU=FLOTM(j,i)*PDTM(j,i);
            %If FLOTM(j,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        if(FSOTM(j+1,i)<=FSCR1)
            FPDYD=PDTM(j+1,i);
        else
            FPDYD=FLOTM(j+1,i)*PDTM(j+1,i);
            %If FLOTM(j,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        
        FPDY(j,i)=-2.0*dtb*(FPDYD-FPDYU)/(dy(j)+dy(j-1));%[Pa.sec/m == kg/m^2/sec]
        
        %solid at cell (j,i), set FPDX=0.0 mandatorily
        if(FLOTM(j,i)<FLErr)
            FPDY(j-1,i)=0.0;%cell up face
            FPDY(j,i)=0.0;%cell bottom face
        end
        
    end
end

%Top & Bottom boundary
for i=1:NIX+2
    if(FSOTM(1,i)<=FSCR1)
        FPDYU=PDTM(1,i);
    else
        FPDYU=FLOTM(1,i)*PDTM(1,i);
        %If FLOTM(1,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSOTM(2,i)<=FSCR1)
        FPDYD=PDTM(2,i);
    else
        FPDYD=FLOTM(2,i)*PDTM(2,i);
        %If FLOTM(2,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    FPDY(1,i)=-2.0*dtb*(FPDYD-FPDYU)/(2.0*dy(1));%upper boundary [Pa.sec/m == kg/m^2/sec]
    %NOTE: FLOTM(1,1:NIX+2)=FLOTM(2,1:NIX+2),
    %PDTM(1,1:NIX+2)=PDTM(2,1:NIX+2) --> FPDY(1,1:NIX+2)=0.0
    
    %solid at cell (2,i), set FPDX=0.0 mandatorily
    if(FLOTM(2,i)<FLErr)
        FPDY(1,i)=0.0;%cell up face
        FPDY(2,i)=0.0;%cell bottom face
    end
    
    if(FSOTM(NIY+1,i)<=FSCR1)
        FPDYU=PDTM(NIY+1,i);
    else
        FPDYU=FLOTM(NIY+1,i)*PDTM(NIY+1,i);
        %If FLOTM(NIY+1,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSOTM(NIY+2,i)<=FSCR1)
        FPDYD=PDTM(NIY+2,i);
    else
        FPDYD=FLOTM(NIY+2,i)*PDTM(NIY+2,i);
        %If FLOTM(NIY+2,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    FPDY(NIY+1,i)=-2.0*dtb*(FPDYD-FPDYU)/(2.0*dy(NIY));%bottom boundary [Pa.sec/m == kg/m^2/sec]
    %NOTE: FLOTM(NIY+1,1:NIX+2)=FLOTM(NIY+2,1:NIX+2),
    %PDTM(NIY+1,1:NIX+2)=PDTM(NIY+2,1:NIX+2) --> FPDY(NIY+1,1:NIX+2)=0.0

    %solid at cell (NIY+1,i), set FPDX=0.0 mandatorily
    if(FLOTM(NIY+1,i)<=FLErr)
        FPDY(NIY,i)=0.0;%cell up face
        FPDY(NIY+1,i)=0.0;%cell bottom face
    end
    
end

FP0YU=0.0;
FP0YD=0.0;
PSFY=zeros(NIY+1,NIX+2);%Initial static pressure * FLN, to extract static pressure from actual pressure
%NOTE: FL is FLN (New step FL)!
for i=1:NIX+2
    for j=2:NIY
        if(FSNTM(j,i)<=FSCR1)
            FP0YU=P0(j,i);
        else
            FP0YU=FLNTM(j,i)*P0(j,i);
            %If FLNTM(j,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        if(FSNTM(j+1,i)<=FSCR1)
            FP0YD=P0(j+1,i);
        else
            FP0YD=FLNTM(j+1,i)*P0(j+1,i);
            %If FLNTM(j,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
        end
        
        PSFY(j,i)=-2.0*dtb*(FP0YD-FP0YU)/(dy(j)+dy(j-1));%[Pa.sec/m == kg/m^2/sec]
        
        %If VYOTM(j,i)=0.0, effective static pressure does not work, which
        %means that not any flux can go through cell face
        if(FLOTM(j+1,i)<FLErr)
            PSFY(j,i)=0.0;
            PSFY(j+1,i)=0.0;
        end
        
    end
end

%Top & Bottom boundary
for i=1:NIX+2
    if(FSNTM(1,i)<=FSCR1)
        FP0YU=P0(1,i);
    else
        FP0YU=FLNTM(1,i)*P0(1,i);
            %If FLNTM(1,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSNTM(2,i)<=FSCR1)
        FP0YD=P0(2,i);
    else
        FP0YD=FLNTM(2,i)*P0(2,i);
            %If FLNTM(2,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    PSFY(1,i)=-2.0*dtb*(FP0YD-FP0YU)/(2.0*dy(1));%upper boundary --> ==-rho*g*dy [Pa.sec/m == kg/m^2/sec]
    
    if(FLOTM(2,i)<FLErr)
        PSFY(1,i)=0.0;
        PSFY(2,i)=0.0;
    end
    
    if(FSNTM(NIY+1,i)<=FSCR1)
        FP0YU=P0(NIY+1,i);
    else
        FP0YU=FLNTM(NIY+1,i)*P0(NIY+1,i);
            %If FLNTM(NIY+1,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    if(FSNTM(NIY+2,i)<=FSCR1)
        FP0YD=P0(NIY+2,i);
    else
        FP0YD=FLNTM(NIY+2,i)*P0(NIY+2,i);
            %If FLNTM(NIY+2,i)==0.0, then it is boundary inside, no dynamic pressure or liquid pressure at all 
    end
    
    PSFY(NIY+1,i)=-2.0*dtb*(FP0YD-FP0YU)/(2.0*dy(NIY));%bottom boundary --> ==-rho*g*dy [Pa.sec/m == kg/m^2/sec]
    if(FLOTM(NIY+1,i)<FLErr)
        PSFY(NIY,i)=0.0;
        PSFY(NIY+1,i)=0.0;
    end
    
end
%----------------------- 22. y- pressure gradient ---------------------------

%---------------------------- 23. y- gravity --------------------------------

FRGOTM=zeros(NIY+1,NIX+2);%G==gravity related terms in RFVYTM
FRGNTM=zeros(NIY+1,NIX+2);
FROLSU=0.0;%FS*RS+FL*RL, O=Old, L=Liquid, S=Solid, U=Up
FRNLSU=0.0;%FS*RS+FL*RL, N=New, L=Liquid, S=Solid, U=Up
FROLSD=0.0;%FS*RS+FL*RL, O=Old, L=Liquid, S=Solid, D=Down
FRNLSD=0.0;%FS*RS+FL*RL, N=New, L=Liquid, S=Solid, D=Down
% %************** Scheme 1: use approximation, see PPT P.45 *****************
% for i=1:NIX+2
%     for j=1:NIY+1
%         FRGOTM(j,i)=0.5*(FLOTM(j,i)+FLOTM(j+1,i))*(0.5*(RLOGTM(j,i)+RLOGTM(j+1,i))-0.5*(RL0(j,i)+RL0(j+1,i)));
%         FRGNTM(j,i)=0.5*(FLNTM(j,i)+FLNTM(j+1,i))*(0.5*(RLNGTM(j,i)+RLNGTM(j+1,i))-0.5*(RL0(j,i)+RL0(j+1,i)));
%
%         %FRGNTM(j,i)=0.5*(FLNTM(j,i)*RLNGTM(j,i)+FLNTM(j+1,i)*RLNGTM(j+1,i))-0.5*(FLNTM(j,i)*RL0(j,i)+FLNTM(j+1,i)*RL0(j+1,i));
%         %FRGOTM(j,i)=0.5*(FLOTM(j,i)*RLOGTM(j,i)+FLOTM(j+1,i)*RLOGTM(j+1,i))-0.5*(FLOTM(j,i)*RL0(j,i)+FLOTM(j+1,i)*RL0(j+1,i));
%
%         %For VYOYM(j,i)=0.0, this face acts as physical boundary,
%         %RFVYTM(j,i) and corresponding DFLP are set to 0.0 mandatorily
%         if(FLOTM(j+1,i)<FLErr)
%             FRGOTM(j,i)=0.0;
%             FRGOTM(j+1,i)=0.0;
%             FRGNTM(j,i)=0.0;
%             FRGNTM(j+1,i)=0.0;
%         end

%     end
% end
%
% %FRLNTM(j,i)=FRLOTM(j,i);%For constant density. If this is set,
% %then sum(B)=0.0 like lid driven cavity flow
%
% FRG=zeros(NIY+1,NIX+2);
% for i=1:NIX+2
%     for j=1:NIY+1
%         FRG(j,i)=g*dtb*(0.5*FRGOTM(j,i)+0.5*FRGNTM(j,i));%[kg/m^2/sec]
%     end
% end
% %**************************************************************************

%*********** Scheme 2: Use no approximation, see PPT P.44 & P.66 ************
for i=1:NIX+2
    for j=1:NIY+1
        if(FSOTM(j,i)<=FSCR1-1.0e-16)%-1.0e-16 is set just for FSCR1=0.0 so that solid gravity does not show its effect
            FROLSU=FLOTM(j,i)*RLOGTM(j,i)+FSO.OL(j,i)*RSO.OL(j,i)+FSO.OPX(j,i)*RSO.OPX(j,i)+FSO.CPX(j,i)*RSO.CPX(j,i)+FSO.PL(j,i)*RSO.PL(j,i)+FSO.ILM(j,i)*RSO.ILM(j,i);%Up cell with free-moving crystal and liquid at old time step
        else
            FROLSU=FLOTM(j,i)*RLOGTM(j,i);%Up porus cell at old time step
        end
        if(FSOTM(j+1,i)<=FSCR1-1.0e-16)%-1.0e-16 is set just for FSCR1=0.0 so that solid gravity does not show its effect
            FROLSD=FLOTM(j+1,i)*RLOGTM(j+1,i)+FSO.OL(j+1,i)*RSO.OL(j+1,i)+FSO.OPX(j+1,i)*RSO.OPX(j+1,i)+FSO.CPX(j+1,i)*RSO.CPX(j+1,i)+FSO.PL(j+1,i)*RSO.PL(j+1,i)+FSO.ILM(j+1,i)*RSO.ILM(j+1,i);%Down cell with free-moving crystal and liquid at old time step
        else
            FROLSD=FLOTM(j+1,i)*RLOGTM(j+1,i);%Down porus cell at old time step
        end
        if(FSNTM(j,i)<=FSCR1-1.0e-16)%-1.0e-16 is set just for FSCR1=0.0 so that solid gravity does not show its effect
            FRNLSU=FLNTM(j,i)*RLNGTM(j,i)+FS.OL(j,i)*RS.OL(j,i)+FS.OPX(j,i)*RS.OPX(j,i)+FS.CPX(j,i)*RS.CPX(j,i)+FS.PL(j,i)*RS.PL(j,i)+FS.ILM(j,i)*RS.ILM(j,i);%Up cell with free-moving crystal and liquid at new time step
        else
            FRNLSU=FLNTM(j,i)*RLNGTM(j,i);%Up porus cell at new time step
        end
        if(FSNTM(j+1,i)<=FSCR1-1.0e-16)%-1.0e-16 is set just for FSCR1=0.0 so that solid gravity does not show its effect
            FRNLSD=FLNTM(j+1,i)*RLNGTM(j+1,i)+FS.OL(j+1,i)*RS.OL(j+1,i)+FS.OPX(j+1,i)*RS.OPX(j+1,i)+FS.CPX(j+1,i)*RS.CPX(j+1,i)+FS.PL(j+1,i)*RS.PL(j+1,i)+FS.ILM(j+1,i)*RS.ILM(j+1,i);%Down cell with free-moving crystal and liquid at new time step
        else
            FRNLSD=FLNTM(j+1,i)*RLNGTM(j+1,i);%Down porus cell at new time step
        end
        
        FRGOTM(j,i)=0.5*(FROLSU+FROLSD);
        FRGNTM(j,i)=0.5*(FRNLSU+FRNLSD);
        
        %         FRGOTM(j,i)=0.5*(FLOTM(j,i)+FLOTM(j+1,i))*0.5*(RLOGTM(j,i)+RLOGTM(j+1,i));
        %         FRGNTM(j,i)=0.5*(FLNTM(j,i)+FLNTM(j+1,i))*0.5*(RLNGTM(j,i)+RLNGTM(j+1,i));
        
        %For VYOYM(j,i)=0.0, this face acts as physical boundary,
        %RFVYTM(j,i) and corresponding DFLP are set to 0.0 mandatorily
        if(FLOTM(j+1,i)<FLErr)
            FRGOTM(j,i)=0.0;
            FRGOTM(j+1,i)=0.0;
            FRGNTM(j,i)=0.0;
            FRGNTM(j+1,i)=0.0;
        end
        
    end
end

%FRLNTM(j,i)=FRLOTM(j,i);%For constant density. If this is set,
%then sum(B)=0.0 like lid driven cavity flow

FRG=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        FRG(j,i)=g*dtb*0.5*(FRGOTM(j,i)+FRGNTM(j,i))+PSFY(j,i);%[kg/m^2/sec]
        %NOTE: MATLAB is working on 64 bits precision, giving ~16 significant digits. Even in R2018, MATLAB does not support 128 bits precision!
    end
end
% %************************************************************************

%NOTE: In Boussinesq assumption, only change of density in gravity term
%has effect!
%---------------------------- 23. y- gravity --------------------------------

%======================= APPROXIMATION OF RFVXY ============================
%------------------------------ RFVXTM -------------------------------------
RFVXTM=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    for j=2:NIY+1
        RFVXTM(j,i)=(RFVX(j,i)+VXCONV(j-1,i)+VSXDIF(j-1,i)+FPDX(j,i)+PSFX(j,i))/(1.0+dtb*UFRKX(j,i));%+DRFVSXT(j,i)+VSXCONVT(j-1,i) [kg/m^2/sec]
    end
end
RFVXTM(:,1)=0.0;%Proved to be 0.0, see slides p.25-26 in CrystalMove.pptx
RFVXTM(:,NIX+1)=0.0;%Proved to be 0.0, see slides p.25-26 in CrystalMove.pptx
%----------------------------- RFVXTM ------------------------------------

%----------------------------- RFVYTM ------------------------------------
RFVYTM=zeros(NIY+1,NIX+2);
for i=2:NIX+1
    for j=2:NIY
        RFVYTM(j,i)=(RFVY(j,i)+VYCONV(j,i-1)+VSYDIF(j,i-1)+FPDY(j,i)+FRG(j,i))/(1.0+dtb*UFRKY(j,i));%+DRFVSYT(j,i)+VSYCONVT(j,i-1)[kg/m^2/sec]
    end
end

for i=2:NIX+1
    %top impermeable boundary
    RFVYTM(1,i)=(RFVY(1,i)+VYCONV(1,i-1)+VSYDIF(1,i-1)+FPDY(1,i)+FRG(1,i))/(1.0+dtb*UFRKY(1,i));%+DRFVSYT(1,i)+VSYCONVT(1,i-1) This formula gives !=0.0
    
    %bottom impermeable boundary
    RFVYTM(NIY+1,i)=(RFVY(NIY+1,i)+VYCONV(NIY+1,i-1)+VSYDIF(NIY+1,i-1)+FPDY(NIY+1,i)+FRG(NIY+1,i))/(1.0+dtb*UFRKY(NIY+1,i));%+DRFVSYT(NIY+1,i)+VSYCONVT(NIY+1,i-1)This formula gives !=0.0
end
RFVYTM(1,:)=0.0;%Proved to be 0.0, see slides p.25-26 in CrystalMove.pptx
RFVYTM(NIY+1,:)=0.0;%Proved to be 0.0, see slides p.25-26 in CrystalMove.pptx
%--------------------------- RFVYTM -----------------------------------
%===================== APPROXIMATION OF RFVXY =========================

%% ============== MASS CONSERVATION RESIDUAL  ================
%----------------------------- RSCT ------------------------------------
%Solid Density Change Total in MASS CONSERVATION EQUATION
RSCT=zeros(NIY,NIX);%[kg/m^3]
% for i=1:NIX
%     for j=1:NIY
%         RSCT(j,i)=0.5*(RSO.OL(j+1,i+1)+RS.OL(j+1,i+1))*dFS.OL(j+1,i+1)+...
%             0.5*(RSO.OPX(j+1,i+1)+RS.OPX(j+1,i+1))*dFS.OPX(j+1,i+1)+...
%             0.5*(RSO.CPX(j+1,i+1)+RS.CPX(j+1,i+1))*dFS.CPX(j+1,i+1)+...
%             0.5*(RSO.PL(j+1,i+1)+RS.PL(j+1,i+1))*dFS.PL(j+1,i+1)+...
%             0.5*(RSO.ILM(j+1,i+1)+RS.ILM(j+1,i+1))*dFS.ILM(j+1,i+1);
%     end
% end
%NOTE(1): This RSCT is absolutely universal to all conditions, see RSO and RS
%NOTE(2): dFS=dFSSM+dFSLH, i.e., total solid volume fraction;

for i=1:NIX
    for j=1:NIY
        RSCT(j,i)=0.5*(RLOTM(j+1,i+1)+RLNTM(j+1,i+1))*dFS.OL(j+1,i+1)+...
            0.5*(RLOTM(j+1,i+1)+RLNTM(j+1,i+1))*dFS.OPX(j+1,i+1)+...
            0.5*(RLOTM(j+1,i+1)+RLNTM(j+1,i+1))*dFS.CPX(j+1,i+1)+...
            0.5*(RLOTM(j+1,i+1)+RLNTM(j+1,i+1))*dFS.PL(j+1,i+1)+...
            0.5*(RLOTM(j+1,i+1)+RLNTM(j+1,i+1))*dFS.ILM(j+1,i+1);
    end
end
%NOTE(1): This RSCT is Bounssinesq approximation. RSO.ANY and RS.ANY have been set equal to RLO and RLN, while RLNTM, RLOTM have been set as RLB (initial liquid
%density for Boussinesq approximation). Obviously, RSCT + RLNTM*FLNTM-RLOTM*FLOTM should be 0.0!
%----------------------------- RSCT ------------------------------------

%----------------------- x- solid mass flow ----------------------------

%----------------------- x- solid mass flow ----------------------------

%----------------------- y- solid mass flow ----------------------------

%----------------------- y- solid mass flow ----------------------------
 
%------------------------------ MSM ------------------------------------

% %Total Mass of Solid Move in x-axis
 MSMXT=zeros(NIY,NIX);%[kg/m^2/sec]

% %Total Mass of Solid Move in y-axis
 MSMYT=zeros(NIY,NIX);%[kg/m^2/sec]

%------------------------------ MSM ------------------------------------

%------------------------------ RES ------------------------------------
RES=zeros(NIY,NIX);
for i=1:NIX
    for j=1:NIY
        %RES(j,i)=dy(j)*(RFVXTM(j+1,i+1)-RFVXTM(j+1,i))+dx(i)*(RFVYTM(j+1,i+1)-RFVYTM(j,i+1))+RSCT(j,i)*dx(i)*dy(j)/dtb+(FLNTM(j+1,i+1)*RLNTM(j+1,i+1)-FLOTM(j+1,i+1)*RLOTM(j+1,i+1))*dx(i)*dy(j)/dtb+dy(j)*MSMXT(j,i)+dx(i)*MSMYT(j,i);%[kg/m/sec == Pa.sec]
        %NOTE: This RES is absolutely universal to all conditions, especially for those related to significant density change. However, in this
        %function, input parameters RLNTM, RLOTM, RSO.ANY and RS.ANY have been set as RLB (initial liquid density for Boussinesq approximation). In VPGM.m function,
        %RLNTM, RLOTM, RSO.ANY and RS.ANY are n+1 and n step liquid and solid density, respectively ,and may be different from each other.
        
        RES(j,i)=dy(j)*(RFVXTM(j+1,i+1)-RFVXTM(j+1,i))+dx(i)*(RFVYTM(j+1,i+1)-RFVYTM(j,i+1))+RSCT(j,i)*dx(i)*dy(j)/dtb+(FLNTM(j+1,i+1)*RLNTM(j+1,i+1)-FLOTM(j+1,i+1)*RLOTM(j+1,i+1))*dx(i)*dy(j)/dtb+dy(j)*MSMXT(j,i)+dx(i)*MSMYT(j,i);%[kg/m/sec == Pa.sec]
        %NOTE: This RES is Bounssinesq approximation. RLOTM, RLNTM, RSO.ANY, RS.ANY are RLB.
        
    end
end
%Detailed information can be accessed in PPT.

for i=1:NIX
    B(1+(i-1)*NIY:i*NIY)=-RES(1:NIY,i)/dtb;%[kg/m/sec == Pa.sec]
end
%B(1)=0.0;
B(NIX*NIY)=0.0;

% RHS=sum(B);
% [u,s]=eig(A);%s is eigenvalue which is used to determine wether matrix A is positive-definite or not
%if A(1,1) is not set as 1.0 etc., s has both positive and negative values, though negative
%value is rather small which is can be taken as 0. Thus, matrix A here
%is positive-semidefinite.

%------------------------------ RES ------------------------------------
%=================== MASS CONSERVATION RESIDUAL  ======================

%% ================= PRESSURE CORRECTION =====================
%Solve linear equations [A][X]=[B]
DFLP=A\B;%[kg/m/sec^2 == Pa]

% AS=sparse(A);%convert sparse A into sparse storage format
% DFLP=PardisoSolver(AS,B);
%NOTE: Any warning, error message from PardisoSolver can be looked up in
%Pardiso_manual. Mostly, error comes from liscence .lic to whom .lic is
%valid.

%     DFLP=GaussSeidel(A,B);
%     DFLP2=pentaDiag_solve(A,B);
%     DFLP=DFLP2-DFLP1;
%     [DFLP,errr,iterr,flag]=sor(A,B,B,0.5,2000,1e-4);
%     DFLP=DFLP1;
%     DFLP=conjgrad(A,B);
%     DFLP=ConjuGrad(A,B,POTM-P0);
%     [DFLP,flags]=PCG(A,B);

%---------------------- New  Relative Pressure ------------------------
%Interior grids
for i=2:NIX+1
    for j=2:NIY+1
        k=(i-2)*NIY+j-1;
        if(FSNTM(j,i)<=FSCR1)%Crystal Free-Moving Region
            PDTM(j,i)=PDTM(j,i)+DFLP(k);
        else%Porus Region, Crystals are interconnected
            PDTM(j,i)=(PDTM(j,i)*FLOTM(j,i)+DFLP(k))/FLNTM(j,i);
            if(FLNTM(j,i)<=FLErr)%this cell is complete solid, we set it mandatorily
                PDTM(j,i)=1.0;
            end
        end
    end
end

for i=2:NIY+1
    %left boundary: VX(1:NIY+2,1)=0, RFVXTM(1:NIY+2,1)=0, so DFLP(1:NIY)=DFLP_left_ghost_cell
    if(FSNTM(i,1)<=FSCR1)%Crystal Free-Moving Region
        PDTM(i,1)=PDTM(i,1)+DFLP(i-1);
    else%Porus Region, Crystals are interconnected
        
        PDTM(i,1)=(PDTM(i,1)*FLOTM(i,1)+DFLP(i-1))/FLNTM(i,1);
        %solid cell, set it mandatorily
        if(FLNTM(i,1)<=FLErr)
            PDTM(i,1)=1.0;
        end
    end
    
    %right boundary: VX(1:NIY+2,NIX+1)=0, RFVXTM(1:NIY+2,NIX+1)=0, so DFLP((NIX-1)*NIY+1:NIX*NIY)=DFLP_righ_ghost_cell
    if(FSNTM(i,NIX+2)<=FSCR1)%Crystal Free-Moving Region
        PDTM(i,NIX+2)=PDTM(i,NIX+2)+DFLP((NIX-1)*NIY+i-1);
    else%Porus Region, Crystals are interconnected
        
        PDTM(i,NIX+2)=(PDTM(i,NIX+2)*FLOTM(i,NIX+2)+DFLP((NIX-1)*NIY+i-1))/FLNTM(i,NIX+2);
        %solid cell, set it mandatorily
        if(FLNTM(i,NIX+2)<=FLErr)
            PDTM(i,NIX+2)=1.0;
        end
    end
end

for i=1:NIX
    %top boundary: DFLP(1:NIY:(NIX-1)*NIY+1)=DFLP_upper_ghost_cell
    %VERY IMPORTANT NOTE: IT IS NOT PDTM(1,i+1)-DFLP(1+(i-1)*NIY), or PDTM(1,i+1)*FLOTM(1,i+1)-DFLP(1+(i-1)*NIY) because the only driving force in VY (see
    %RFVYTM) is FRG, then in next call of this function FPDY(1,i+1) must be 0.0, meaning PATM(2,i+1)_last_call==PDTM(1,i+1)_last_call. So, we have to set 
    %PDTM(2,i+1)=PDTM(1,i+1) in this call. See also RayleighBenardConvection.m
    if(FSNTM(1,i+1)<=FSCR1)%Crystal Free-Moving Region
        PDTM(1,i+1)=PDTM(1,i+1)+DFLP(1+(i-1)*NIY);
    else%Porus Region, Crystals are interconnected
        
        PDTM(1,i+1)=(PDTM(1,i+1)*FLOTM(1,i+1)+DFLP(1+(i-1)*NIY))/FLNTM(1,i+1);
        %solid cell, set it mandatorily
        if(FLNTM(1,i+1)<=FLErr)
            PDTM(1,i+1)=1.0;
        end
    end
    
    %bottom boundary: DFLP(NIY:NIY:NIX*NIY)=DFLP_bottom_ghost_cell
    %VERY IMPORTANT NOTE: IT IS NOT PDTM(NIY+2,i+1)-DFLP(i*NIY), or PDTM(NIY+2,i+1)*FLOTM(NIY+2,i+1)+DFLP(i*NIY) because the only driving force in VY (see
    %RFVYTM) is FRG, then in next call of this function FPDY(NIY+1,i+1) must be 0.0, meaning PATM(NIY+2,i+1)_last_call==PDTM(NIY+1,i+1)_last_call. So, we have 
    %to set PDTM(NIY+2,i+1)=PDTM(NIY+1,i+1) in this call. See also RayleighBenardConvection.m
    if(FSNTM(NIY+2,i+1)<=FSCR1)%Crystal Free-Moving Region
        PDTM(NIY+2,i+1)=PDTM(NIY+2,i+1)+DFLP(i*NIY);
    else%Porus Region, Crystals are interconnected
        
        PDTM(NIY+2,i+1)=(PDTM(NIY+2,i+1)*FLOTM(NIY+2,i+1)+DFLP(i*NIY))/FLNTM(NIY+2,i+1);
        %solid cell, set it mandatorily
        if(FLNTM(NIY+2,i+1)<=FLErr)
            PDTM(NIY+2,i+1)=1.0;
        end
    end
end
PDTM(1,1)=PDTM(1,2);
PDTM(NIY+2,1)=PDTM(NIY+2,2);
PDTM(1,NIX+2)=PDTM(1,NIX+1);
PDTM(NIY+2,NIX+2)=PDTM(NIY+2,NIX+1);
%---------------------- New Relative Pressure -------------------------

%======================= PRESSURE CORRECTION ==========================

%% =============== MASS FLOW CORRECTION ======================
%-------------------------- DRFVX RFVXN--------------------------------
RFVXN=zeros(NIY+2,NIX+1);%new RFVX
for i=2:NIX
    for j=2:NIY+1
        %VX(1:NIY+2,1)=0.0, RFVXTM(1:NIY+2,1)=0.0 --> DRFVX(1:NIY+2,1)=0.0
        %VX(1:NIY+2,NIX+1)=0.0, RFVXTM(1:NIY+2,NIX+1)=0.0 --> DRFVX(1:NIY+2,NIX+1)=0.0
        k=(i-2)*NIY+j-1;
        DRFVX(j,i)=-2.0*dtb*(DFLP(k+NIY)-DFLP(k))/((dx(i)+dx(i-1))*(1.0+dtb*UFRKX(j,i)));%[Pa.sec/m == kg/m^2/sec]
        RFVXN(j,i)=RFVXTM(j,i)+DRFVX(j,i);
    end
end
%Top NO SLIP boundary
RFVXN(1,1:NIX+1)=-RFVXN(2,1:NIX+1);
%Bottom FREE boundary
RFVXN(NIY+2,:)=RFVXN(NIY+1,:);

% %Top FREE boundary
% RFVXN(1,1:NIX+1)=RFVXN(2,1:NIX+1);
% %Bottom FREE boundary
% RFVXN(NIY+2,:)=RFVXN(NIY+1,:);

%RFVXN(1:NIY+2,1)=0.0
%RFVXN(1:NIY+2,NIX+1)=0.0

%----------------------------- DRFVX ----------------------------------

%-------------------------- DRFVY RFVYN -------------------------------
RFVYN=zeros(NIY+1,NIX+2);%new RFVY
for i=2:NIX+1
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0, RFVYTM(1,1:NIX+2)=0.0 --> DRFVY(1,1:NIX+2)=0.0
        %VY(NIY+1,1:NIX+2)=0.0, RFVYTM(NIY+1,1:NIX+2)=0.0 --> DRFVY(NIY+1,1:NIX+2)=0.0
        k=(i-2)*NIY+j-1;
        DRFVY(j,i)=-2.0*dtb*(DFLP(k+1)-DFLP(k))/((dy(j)+dy(j-1))*(1.0+dtb*UFRKY(j,i)));%[Pa.sec/m == kg/m^2/sec]
        RFVYN(j,i)=RFVYTM(j,i)+DRFVY(j,i);
    end
end
%Left FREE boundary
RFVYN(1:NIY+1,1)=RFVYN(1:NIY+1,2);
%Right FREE bounday
RFVYN(1:NIY+1,NIX+2)=RFVYN(1:NIY+1,NIX+1);

% %Left NO SLIP boundary
% RFVYN(1:NIY+1,1)=-RFVYN(1:NIY+1,2);
% %Right NO SLIP bounday
% RFVYN(1:NIY+1,NIX+2)=-RFVYN(1:NIY+1,NIX+1);

%RFVYN(1,1:NIX+2)=0.0
%RFVYN(NIX+1,1:NIX+2)=0.0

%----------------------------- DRFVY ----------------------------------
%======================= MASS FLOW CORRECTION =========================

%------------------------------ NEW VX ------------------------------------
for i=2:NIX
    for j=2:NIY+1
        %VX(1:NIY+2,1)=0.0
        %VX(1:NIY+2,NIX+1)=0.0
        VXOTM(j,i)=RFVXN(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1)));%[m/sec]
        
        %solid cell, set 0.0 mandatorily
        if(FLNTM(j,i)<=FLErr)
            VXOTM(j,i-1)=0.0;%cell left face
            VXOTM(j,i)=0.0;%cell right face
        end
    end
end

for i=1:NIX+1
    %top NO SLIP boundary
    VXOTM(1,i)=-VXOTM(2,i);
    
    %bottom FREE boundary
    VXOTM(NIY+2,i)=VXOTM(NIY+1,i);
end

% for i=1:NIX+1
%     %top FREE boundary
%     VXOTM(1,i)=VXOTM(2,i);
%     
%     %bottom FREE boundary
%     VXOTM(NIY+2,i)=VXOTM(NIY+1,i);
% end

%------------------------------ NEW VX ------------------------------------

%------------------------------ NEW VY ------------------------------------
for i=2:NIX+1
    %VY(1,1:NIX+2)=0.0 
    %VY(NIY+1,1:NIX+2)=0.0 
    for j=2:NIY
        VYOTM(j,i)=RFVYN(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i)));
        
        %solid cell, set 0.0 mandatorily
        if(FLNTM(j,i)<=FLErr)
            VYOTM(j-1,i)=0.0;%cell up face
            VYOTM(j,i)=0.0;%cell bottom face
        end
        
    end
end

for i=1:NIY+1
    %left FREE boundary
    VYOTM(i,1)=VYOTM(i,2);
    
    %right FREE boundary
    VYOTM(i,NIX+2)=VYOTM(i,NIX+1);
end

% for i=1:NIY+1
%     %left NO SLIP boundary
%     VYOTM(i,1)=-VYOTM(i,2);
%
%     %right NO SLIP boundary
%     VYOTM(i,NIX+2)=-VYOTM(i,NIX+1);
% end
%------------------------------ NEW VY ------------------------------------

%% ===================== MASS BALANCE CHECK ======================
RESM=zeros(NIY+2,NIX+2);%absolute residual of mass
RESMR=zeros(NIY+2,NIX+2);%relative residual of mass
mass_check=1.0;%mass check marker
for i=2:NIX+1
    for j=2:NIY+1
        % RESM(j,i)=dtb*dy(j-1)*(RFVXN(j,i)-RFVXN(j,i-1))+dtb*dx(i-1)*(RFVYN(j,i)-RFVYN(j-1,i))+RSCT(j-1,i-1)*dx(i-1)*dy(j-1)+(FLNTM(j,i)*RLNTM(j,i)-FLOTM(j,i)*RLOTM(j,i))*dx(i-1)*dy(j-1)+dtb*dy(j-1)*MSMXT(j-1,i-1)+dtb*dx(i-1)*MSMYT(j-1,i-1);%absolute mass residual [kg/sec]
        %NOTE: This RESM is absolutely universal to all conditions, especially for those related to significant density change. However, in this
        %function, input parameters RLNTM, RLOTM, RSO.ANY and RS.ANY have been set as RLB (initial liquid density for Boussinesq approximation). In VPGM.m function,
        %RLNTM, RLOTM, RSO.ANY and RS.ANY are n+1 and n step liquid and solid density, respectively ,and may be different from each other.
        
        RESM(j,i)=dtb*dy(j-1)*(RFVXN(j,i)-RFVXN(j,i-1))+dtb*dx(i-1)*(RFVYN(j,i)-RFVYN(j-1,i))+RSCT(j-1,i-1)*dx(i-1)*dy(j-1)+(FLNTM(j,i)*RLNTM(j,i)-FLOTM(j,i)*RLOTM(j,i))*dx(i-1)*dy(j-1)+dtb*dy(j-1)*MSMXT(j-1,i-1)+dtb*dx(i-1)*MSMYT(j-1,i-1);%absolute mass residual [kg/sec]
        %NOTE: This RES is Bounssinesq approximation. RLOTM, RLNTM, RSO.ANY, RS.ANY are RLB.
        
        RESMR(j,i)=RESM(j,i)/(dx(i-1)*dy(j-1)*RL0(j,i)/dtb);%relative mass residual with respect to initial mass flow [(kg/sec)/(kg/sec)==1]
        if(abs(RESM(j,i))>1.0e-9)
            fprintf(2,'Mass check: (%2d,%2d) not balanced! %E\n',j-1,i-1,RESM(j,i));
            mass_check=-1.0;
        end
        
    end
end

RESM(1,2:NIX+1)=RESM(2,2:NIX+1);
RESM(NIY+2,2:NIX+1)=RESM(NIY+1,2:NIX+1);
RESM(1:NIY+2,1)=RESM(1:NIY+2,2);
RESM(1:NIY+2,NIX+2)=RESM(1:NIY+2,NIX+1);

if(mass_check>0.0)
    fprintf('Mass Balanced -- Max error %E\n',max(max(abs(RESM))));
end
%------------------------- MASS BALANCE CHECK -----------------------------

%% ================== RETURN NEW VARIABLES ======================
P=PDTM;
VX=VXOTM;
VY=VYOTM;
%------------------------ RETURN NEW VARIABLES ----------------------------

end
