function NASV(VXTM,VYTM,MUNTM,RLNTM)
%This function is used to to calculate n+1 step absolute solid velocity, from Eq.(2-13 ~ 2-15) Wangshenqiang
%Created on 2020-7-21

%========== INPUTS ==========
%VXTM: old absolute liquid x-axis velocity [m/sec]
%VYTM: old absolute liquid y-axis velocity [m/sec]
%MUNTM: new viscosity dynamic of pure liquid [Pa.sec]
%RLNTM: new liquid density [kg/m^3]
%VSXR: old relative solid x-axis velocity [m/sec]
%VSYR: old relative solid y-axis velocity [m/sec]
%RS: new solid density [kg/m^3]
%FS: new solid volume fraction [1]

%========= OUTPUTS ==========
%VSX.ANY: new absolute solid x-axis velosity [m/sec]
%VSY.ANY: new absolute solid y-axis velosity [m/sec]

global NIX
global NIY
global dx
global dy
global dtb
global VSXR
global VSYR
global RS
global FS
global dFSLH
global g
global FSCR1
global VSX
global VSY

%% =================== INITIALIZE VSXR VSYR ==================
%Any solid crystallized due to cooling are assumed 0.0 velocity with respect to x-y-z coordinate system, thus this part solid has a -VX and -VY
%relative velocity!
for i=1:NIX+1
    for j=1:NIY+2
        VSXR(j,i)=-VXTM(j,i);
    end
end

for i=1:NIX+2
    for j=1:NIY+1
        VSYR(j,i)=-VYTM(j,i);
    end
end

%% =================== UPDATE SOLID RADIUS ===================
RP=struct('OL',zeros(NIY+2,NIX+2),...%olivine radius [m]
    'OPX',zeros(NIY+2,NIX+2),...%opx radius [m]
    'CPX',zeros(NIY+2,NIX+2),...%cpx radius [m]
    'PL',zeros(NIY+2,NIX+2),...%pl radius [m]
    'ILM',zeros(NIY+2,NIX+2));%ilmentite radius [m]

% %3D with dz=unit length
% for i=2:NIX+1
%     for j=2:NIY+1
%         RP.OL(j,i)=nthroot((3.0*dFSLH.OL(j,i)*dx(i-1)*dy(j-1)*1.0)/(4.0*pi),3);%dz=unit length
%         RP.OPX(j,i)=nthroot((3.0*dFSLH.OPX(j,i)*dx(i-1)*dy(j-1)*1.0)/(4.0*pi),3);%dz=unit length
%         RP.CPX(j,i)=nthroot((3.0*dFSLH.CPX(j,i)*dx(i-1)*dy(j-1)*1.0)/(4.0*pi),3);%dz=unit length
%         RP.PL(j,i)=nthroot((3.0*dFSLH.PL(j,i)*dx(i-1)*dy(j-1)*1.0)/(4.0*pi),3);%dz=unit length
%         RP.ILM(j,i)=nthroot((3.0*dFSLH.ILM(j,i)*dx(i-1)*dy(j-1)*1.0)/(4.0*pi),3);%dz=unit length
%     end
% end

%2D, assuming same area
for i=2:NIX+1
    for j=2:NIY+1
        RP.OL(j,i)=sqrt(dFSLH.OL(j,i)*dx(i-1)*dy(j-1)/pi);
        RP.OPX(j,i)=sqrt(dFSLH.OPX(j,i)*dx(i-1)*dy(j-1)/pi);
        RP.CPX(j,i)=sqrt(dFSLH.CPX(j,i)*dx(i-1)*dy(j-1)/pi);
        RP.PL(j,i)=sqrt(dFSLH.PL(j,i)*dx(i-1)*dy(j-1)/pi);
        RP.ILM(j,i)=sqrt(dFSLH.ILM(j,i)*dx(i-1)*dy(j-1)/pi);
    end
end

%Olivine boundary
RP.OL(2:NIY+1,1)=RP.OL(2:NIY+1,2);
RP.OL(2:NIY+1,NIX+2)=RP.OL(2:NIY+1,NIX+1);
RP.OL(1,1:NIX+2)=RP.OL(2,1:NIX+2);
RP.OL(NIY+2,1:NIX+2)=RP.OL(NIY+1,1:NIX+2);

%OPX boundary
RP.OPX(2:NIY+1,1)=RP.OPX(2:NIY+1,2);
RP.OPX(2:NIY+1,NIX+2)=RP.OPX(2:NIY+1,NIX+1);
RP.OPX(1,1:NIX+2)=RP.OPX(2,1:NIX+2);
RP.OPX(NIY+2,1:NIX+2)=RP.OPX(NIY+1,1:NIX+2);

%CPX boundary
RP.CPX(2:NIY+1,1)=RP.CPX(2:NIY+1,2);
RP.CPX(2:NIY+1,NIX+2)=RP.CPX(2:NIY+1,NIX+1);
RP.CPX(1,1:NIX+2)=RP.CPX(2,1:NIX+2);
RP.CPX(NIY+2,1:NIX+2)=RP.CPX(NIY+1,1:NIX+2);

%Plagioclase boundary
RP.PL(2:NIY+1,1)=RP.PL(2:NIY+1,2);
RP.PL(2:NIY+1,NIX+2)=RP.PL(2:NIY+1,NIX+1);
RP.PL(1,1:NIX+2)=RP.PL(2,1:NIX+2);
RP.PL(NIY+2,1:NIX+2)=RP.PL(NIY+1,1:NIX+2);

%Ilmenite boundary
RP.ILM(2:NIY+1,1)=RP.ILM(2:NIY+1,2);
RP.ILM(2:NIY+1,NIX+2)=RP.ILM(2:NIY+1,NIX+1);
RP.ILM(1,1:NIX+2)=RP.ILM(2,1:NIX+2);
RP.ILM(NIY+2,1:NIX+2)=RP.ILM(NIY+1,1:NIX+2);

%% =================== UPDATE VSXR.ANY ====================
%New time step solid relative velocity x-axis: VSXR [m/sec]
for i=2:NIX+1
    for j=2:NIY+1
        if((FS.OL(j,i)>0.0)&&(FS.OL(j,i)<=FSCR1))
        VSXR.OL(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.OL(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.OL(j,i)*(RS.OL(j,i)+0.5*RLNTM(j,i))))*VSXR.OL(j,i);%dz=unit length
        else
            VSXR.OL(j,i)=0.0;
        end
        if((FS.OPX(j,i)>0.0)&&(FS.OPX(j,i)<=FSCR1))
        VSXR.OPX(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.OPX(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.OPX(j,i)*(RS.OPX(j,i)+0.5*RLNTM(j,i))))*VSXR.OPX(j,i);%dz=unit length
        else
            VSXR.OPX(j,i)=0.0;
        end
        if((FS.CPX(j,i)>0.0)&&(FS.CPX(j,i)<=FSCR1))
        VSXR.CPX(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.CPX(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.CPX(j,i)*(RS.CPX(j,i)+0.5*RLNTM(j,i))))*VSXR.CPX(j,i);%dz=unit length
        else
            VSXR.CPX(j,i)=0.0;
        end
        if((FS.PL(j,i)>0.0)&&(FS.PL(j,i)<=FSCR1))
        VSXR.PL(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.PL(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.PL(j,i)*(RS.PL(j,i)+0.5*RLNTM(j,i))))*VSXR.PL(j,i);%dz=unit length
        else
            VSXR.PL(j,i)=0.0;
        end
        if((FS.ILM(j,i)>0.0)&&(FS.ILM(j,i)<=FSCR1))
        VSXR.ILM(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.ILM(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.ILM(j,i)*(RS.ILM(j,i)+0.5*RLNTM(j,i))))*VSXR.ILM(j,i);%dz=unit length
        else
            VSXR.ILM(j,i)=0.0;
        end
        
        %IMPORTANT NOTE(1): the reason that dx, dy, RP.ANY, MUOTM, FSTM.ANY, RS.ANY, RLTM are used as general grid values not velocity grid (i.e., x-axis
        %velocity, y-axis velocity, they are displaced half grid along x-axis and y-axis respectively), is that the deduction of vsxr, vsyr and
        %vszr is carried out within one certain cell, that is, all information used to calculate vsxr, vsyr and vszr are that of one certain cell,
        %thus these information are the same. Here, we set dx, dy, dz, RP,ANY, MUOTM, FSTM.ANY, RLTM as general grid values. The relation between
        %vsxr, vsyr, vszr and vsx, vsy, vsz are then described by Eq.(2-20 ~ 2-25) in Wangshenqiang.
        
        %IMPORTANT NOTE(2): all variables are new time step except VSXR, VSYR!
    end
end

%Left impermeable boundary
VSXR.OL(1:NIY+2,1)=-VSXR.OL(1:NIY+2,2);
VSXR.OPX(1:NIY+2,1)=-VSXR.OPX(1:NIY+2,2);
VSXR.CPX(1:NIY+2,1)=-VSXR.CPX(1:NIY+2,2);
VSXR.PL(1:NIY+2,1)=-VSXR.PL(1:NIY+2,2);
VSXR.ILM(1:NIY+2,1)=-VSXR.ILM(1:NIY+2,2);

%Right impermeable boundary
VSXR.OL(1:NIY+2,NIX+2)=-VSXR.OL(1:NIY+2,NIX+1);
VSXR.OPX(1:NIY+2,NIX+2)=-VSXR.OPX(1:NIY+2,NIX+1);
VSXR.CPX(1:NIY+2,NIX+2)=-VSXR.CPX(1:NIY+2,NIX+1);
VSXR.PL(1:NIY+2,NIX+2)=-VSXR.PL(1:NIY+2,NIX+1);
VSXR.ILM(1:NIY+2,NIX+2)=-VSXR.ILM(1:NIY+2,NIX+1);

%Top NO SLIP boundary
VSXR.OL(1,2:NIX+1)=-VSXR.OL(2,2:NIX+1);
VSXR.OPX(1,2:NIX+1)=-VSXR.OPX(2,2:NIX+1);
VSXR.CPX(1,2:NIX+1)=-VSXR.CPX(2,2:NIX+1);
VSXR.PL(1,2:NIX+1)=-VSXR.PL(2,2:NIX+1);
VSXR.ILM(1,2:NIX+1)=-VSXR.ILM(2,2:NIX+1);

% %Top FREE boundary
% VSXR.OL(1,2:NIX+1)=VSXR.OL(2,2:NIX+1);
% VSXR.OPX(1,2:NIX+1)=VSXR.OPX(2,2:NIX+1);
% VSXR.CPX(1,2:NIX+1)=VSXR.CPX(2,2:NIX+1);
% VSXR.PL(1,2:NIX+1)=VSXR.PL(2,2:NIX+1);
% VSXR.ILM(1,2:NIX+1)=VSXR.ILM(2,2:NIX+1);

%Bottom NO SLIP boundary
VSXR.OL(NIY+2,2:NIX+1)=-VSXR.OL(NIY+1,2:NIX+1);
VSXR.OPX(NIY+2,2:NIX+1)=-VSXR.OPX(NIY+1,2:NIX+1);
VSXR.CPX(NIY+2,2:NIX+1)=-VSXR.CPX(NIY+1,2:NIX+1);
VSXR.PL(NIY+2,2:NIX+1)=-VSXR.PL(NIY+1,2:NIX+1);
VSXR.ILM(NIY+2,2:NIX+1)=-VSXR.ILM(NIY+1,2:NIX+1);

% %Bottom FREE boundary
% VSXR.OL(NIY+2,2:NIX+1)=VSXR.OL(NIY+1,2:NIX+1);
% VSXR.OPX(NIY+2,2:NIX+1)=VSXR.OPX(NIY+1,2:NIX+1);
% VSXR.CPX(NIY+2,2:NIX+1)=VSXR.CPX(NIY+1,2:NIX+1);
% VSXR.PL(NIY+2,2:NIX+1)=VSXR.PL(NIY+1,2:NIX+1);
% VSXR.ILM(NIY+2,2:NIX+1)=VSXR.ILM(NIY+1,2:NIX+1);

%% =================== UPDATE VSYR.ANY ====================
%New time step solid relative velocity y-axis: VSYR
for i=2:NIX+1
    for j=2:NIY+1
        if((FS.OL(j,i)>0.0)&&(FS.OL(j,i)<=FSCR1))
        VSYR.OL(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.OL(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.OL(j,i)*(RS.OL(j,i)+0.5*RLNTM(j,i))))*VSYR.OL(j,i)+(RS.OL(j,i)-RLNTM(j,i))*dtb*g/(RS.OL(j,i)+0.5*RLNTM(j,i));%dz=unit length
        else
            VSYR.OL(j,i)=0.0;
        end
        if((FS.OPX(j,i)>0.0)&&(FS.OPX(j,i)<=FSCR1))
        VSYR.OPX(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.OPX(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.OPX(j,i)*(RS.OPX(j,i)+0.5*RLNTM(j,i))))*VSYR.OPX(j,i)+(RS.OPX(j,i)-RLNTM(j,i))*dtb*g/(RS.OPX(j,i)+0.5*RLNTM(j,i));
        else
            VSYR.OPX(j,i)=0.0;
        end
        if((FS.CPX(j,i)>0.0)&&(FS.CPX(j,i)<=FSCR1))
        VSYR.CPX(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.CPX(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.CPX(j,i)*(RS.CPX(j,i)+0.5*RLNTM(j,i))))*VSYR.CPX(j,i)+(RS.CPX(j,i)-RLNTM(j,i))*dtb*g/(RS.CPX(j,i)+0.5*RLNTM(j,i));
        else
            VSYR.CPX(j,i)=0.0;
        end
        if((FS.PL(j,i)>0.0)&&(FS.PL(j,i)<=FSCR1))
        VSYR.PL(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.PL(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.PL(j,i)*(RS.PL(j,i)+0.5*RLNTM(j,i))))*VSYR.PL(j,i)+(RS.PL(j,i)-RLNTM(j,i))*dtb*g/(RS.PL(j,i)+0.5*RLNTM(j,i));
        else
            VSYR.PL(j,i)=0.0;
        end
        if((FS.ILM(j,i)>0.0)&&(FS.ILM(j,i)<=FSCR1))
        VSYR.ILM(j,i)=(1.0-6.0*pi*MUNTM(j,i)*RP.ILM(j,i)*dtb/(dx(i-1)*dy(j-1)*1.0*FS.ILM(j,i)*(RS.ILM(j,i)+0.5*RLNTM(j,i))))*VSYR.ILM(j,i)+(RS.ILM(j,i)-RLNTM(j,i))*dtb*g/(RS.ILM(j,i)+0.5*RLNTM(j,i));
        else
            VSYR.ILM(j,i)=0.0;
        end
        
        %IMPORTANT NOTE(1): the reason that dx, dy, RP.ANY, MUOTM, FSTM.ANY, RS.ANY, RLTM are used as general grid values not velocity grid (i.e., x-axis
        %velocity, y-axis velocity, they are displaced half grid along x-axis and y-axis respectively), is that the deduction of vsxr, vsyr and
        %vszr is carried out within one certain cell, that is, all information used to calculate vsxr, vsyr and vszr are that of one certain cell,
        %thus these information are the same. Here, we set dx, dy, dz, RP,ANY, MUOTM, FSTM.ANY, RLTM as general grid values. The relation between
        %vsxr, vsyr, vszr and vsx, vsy, vsz are then described by Eq.(2-20 ~ 2-25) in Wangshenqiang.
        
        %IMPORTANT NOTE(2): all variables are new time step except VSXR, VSYR!
        
        %IMPORTANT NOTE(3): it's + not - in y-axis since y axis is downward in our case.
    end
end

% %Left NO SLIP boundary
% VSYR.OL(2:NIY+1,1)=-VSYR.OL(2:NIY+1,2);
% VSYR.OPX(2:NIY+1,1)=-VSYR.OPX(2:NIY+1,2);
% VSYR.CPX(2:NIY+1,1)=-VSYR.CPX(2:NIY+1,2);
% VSYR.PL(2:NIY+1,1)=-VSYR.PL(2:NIY+1,2);
% VSYR.ILM(2:NIY+1,1)=-VSYR.ILM(2:NIY+1,2);

%Left FREE boundary
VSYR.OL(2:NIY+1,1)=VSYR.OL(2:NIY+1,2);
VSYR.OPX(2:NIY+1,1)=VSYR.OPX(2:NIY+1,2);
VSYR.CPX(2:NIY+1,1)=VSYR.CPX(2:NIY+1,2);
VSYR.PL(2:NIY+1,1)=VSYR.PL(2:NIY+1,2);
VSYR.ILM(2:NIY+1,1)=VSYR.ILM(2:NIY+1,2);

% %Right NO SLIP boundary
% VSYR.OL(2:NIY+1,NIX+2)=-VSYR.OL(2:NIY+1,NIX+1);
% VSYR.OPX(2:NIY+1,NIX+2)=-VSYR.OPX(2:NIY+1,NIX+1);
% VSYR.CPX(2:NIY+1,NIX+2)=-VSYR.CPX(2:NIY+1,NIX+1);
% VSYR.PL(2:NIY+1,NIX+2)=-VSYR.PL(2:NIY+1,NIX+1);
% VSYR.ILM(2:NIY+1,NIX+2)=-VSYR.ILM(2:NIY+1,NIX+1);

%Right FREE boundary
VSYR.OL(2:NIY+1,NIX+2)=VSYR.OL(2:NIY+1,NIX+1);
VSYR.OPX(2:NIY+1,NIX+2)=VSYR.OPX(2:NIY+1,NIX+1);
VSYR.CPX(2:NIY+1,NIX+2)=VSYR.CPX(2:NIY+1,NIX+1);
VSYR.PL(2:NIY+1,NIX+2)=VSYR.PL(2:NIY+1,NIX+1);
VSYR.ILM(2:NIY+1,NIX+2)=VSYR.ILM(2:NIY+1,NIX+1);

%Top impermeable boundary
VSYR.OL(1,1:NIX+2)=-VSYR.OL(2,1:NIX+2);
VSYR.OPX(1,1:NIX+2)=-VSYR.OPX(2,1:NIX+2);
VSYR.CPX(1,1:NIX+2)=-VSYR.CPX(2,1:NIX+2);
VSYR.PL(1,1:NIX+2)=-VSYR.PL(2,1:NIX+2);
VSYR.ILM(1,1:NIX+2)=-VSYR.ILM(2,1:NIX+2);

%Bottom impermeable boundary
VSYR.OL(NIY+2,1:NIX+2)=-VSYR.OL(NIY+1,1:NIX+2);
VSYR.OPX(NIY+2,1:NIX+2)=-VSYR.OPX(NIY+1,1:NIX+2);
VSYR.CPX(NIY+2,1:NIX+2)=-VSYR.CPX(NIY+1,1:NIX+2);
VSYR.PL(NIY+2,1:NIX+2)=-VSYR.PL(NIY+1,1:NIX+2);
VSYR.ILM(NIY+2,1:NIX+2)=-VSYR.ILM(NIY+1,1:NIX+2);

%% =================== UPDATE VSX.ANY =====================

VSXU=0.0;%Temporary absolute VSX corresponding to x-axis Upwind scheme
VSXD=0.0;%Temporary absolute VSX corresponding to x-axis Downwind scheme

%NOTE: The following will give VSX.ANY=VXTM when FS.ANY=0, this wouldn't
%bother, because VSX.ANY will be multiplied by FS.ANY in all functions, and
%thus solid will not contribute to fluxes.
for i=2:NIX
    for j=2:NIY+1
        %OL x-axis absolute velocity
        VSXU=VXTM(j,i)+VSXR.OL(j,i);
        VSXD=VXTM(j,i)+VSXR.OL(j,i+1);
        if(VSXU>=0.0)
            VSX.OL(j,i)=VSXU;
        end
        if(VSXD<0.0)
            VSX.OL(j,i)=VSXD;
        end
        
        %OPX x-axis absolute velocity
        VSXU=VXTM(j,i)+VSXR.OPX(j,i);
        VSXD=VXTM(j,i)+VSXR.OPX(j,i+1);
        if(VSXU>=0.0)
            VSX.OPX(j,i)=VSXU;
        end
        if(VSXD<0.0)
            VSX.OPX(j,i)=VSXD;
        end
        
        %CPX x-axis absolute velocity
        VSXU=VXTM(j,i)+VSXR.CPX(j,i);
        VSXD=VXTM(j,i)+VSXR.CPX(j,i+1);
        if(VSXU>=0.0)
            VSX.CPX(j,i)=VSXU;
        end
        if(VSXD<0.0)
            VSX.CPX(j,i)=VSXD;
        end
        
        %PL x-axis absolute velocity
        VSXU=VXTM(j,i)+VSXR.PL(j,i);
        VSXD=VXTM(j,i)+VSXR.PL(j,i+1);
        if(VSXU>=0.0)
            VSX.PL(j,i)=VSXU;
        end
        if(VSXD<0.0)
            VSX.PL(j,i)=VSXD;
        end
        
        %ILM x-axis absolute velocity
        VSXU=VXTM(j,i)+VSXR.ILM(j,i);
        VSXD=VXTM(j,i)+VSXR.ILM(j,i+1);
        if(VSXU>=0.0)
            VSX.ILM(j,i)=VSXU;
        end
        if(VSXD<0.0)
            VSX.ILM(j,i)=VSXD;
        end
        
        %NOTE: VSXR.ANY are on general grid!
    end
end

for i=1:NIY+2
    %VSX.ANY(1:NIY+2,1)=0.0 --> left impermeable boundary
    VSX.OL(i,1)=0.0;
    VSX.OPX(i,1)=0.0;
    VSX.CPX(i,1)=0.0;
    VSX.PL(i,1)=0.0;
    VSX.ILM(i,1)=0.0;
    
    %VSX.ANY(1:NIY+2,NIX+1)=0.0 --> right impermeable boundary
    VSX.OL(i,NIX+1)=0.0;
    VSX.OPX(i,NIX+1)=0.0;
    VSX.CPX(i,NIX+1)=0.0;
    VSX.PL(i,NIX+1)=0.0;
    VSX.ILM(i,NIX+1)=0.0;
end

%Top & bottom NO SLIP boundary
for i=1:NIX+1
    %VSX.ANY(1,1:NIX+1)=-VSX.ANY(2,1:NIX+1) --> top NO SLIP boundary
    VSX.OL(1,i)=-VSX.OL(2,i);
    VSX.OPX(1,i)=-VSX.OPX(2,i);
    VSX.CPX(1,i)=-VSX.CPX(2,i);
    VSX.PL(1,i)=-VSX.PL(2,i);
    VSX.ILM(1,i)=-VSX.ILM(2,i);
    
    %VSX.ANY(NIY+2,1:NIX+1)=-VSX.ANY(NIY+1,1:NIX+1) --> bottom NO SLIP boundary
    VSX.OL(NIY+2,i)=-VSX.OL(NIY+1,i);
    VSX.OPX(NIY+2,i)=-VSX.OPX(NIY+1,i);
    VSX.CPX(NIY+2,i)=-VSX.CPX(NIY+1,i);
    VSX.PL(NIY+2,i)=-VSX.PL(NIY+1,i);
    VSX.ILM(NIY+2,i)=-VSX.ILM(NIY+1,i);
end

% Top & Bottom FREE boundary
% for i=1:NIX+1
%     %VSX.ANY(1,1:NIX+1)=VSX.ANY(2,1:NIX+1) --> top FREE boundary
%     VSX.OL(1,i)=VSX.OL(2,i);
%     VSX.OPX(1,i)=VSX.OPX(2,i);
%     VSX.CPX(1,i)=VSX.CPX(2,i);
%     VSX.PL(1,i)=VSX.PL(2,i);
%     VSX.ILM(1,i)=VSX.ILM(2,i);
%     
%     %VSX.ANY(NIY+2,1:NIX+1)=VSX.ANY(NIY+1,1:NIX+1) --> bottom FREE boundary
%     VSX.OL(NIY+2,i)=VSX.OL(NIY+1,i);
%     VSX.OPX(NIY+2,i)=VSX.OPX(NIY+1,i);
%     VSX.CPX(NIY+2,i)=VSX.CPX(NIY+1,i);
%     VSX.PL(NIY+2,i)=VSX.PL(NIY+1,i);
%     VSX.ILM(NIY+2,i)=VSX.ILM(NIY+1,i);
% end

%% =================== UPDATE VSY.ANY =====================

VSYU=0.0;%Temporary absolute VSY correspoding to y-axis Upwind scheme
VSYD=0.0;%Temporary absolute VSY correspoding to y-axis Downwind scheme

%NOTE: The following will give VSX.ANY=VXTM when FS.ANY=0, this wouldn't
%bother, because VSX.ANY will be multiplied by FS.ANY in all functions, and
%thus solid will not contribute to fluxes.
for i=2:NIX+1
    for j=2:NIY
        %OL y-axis absolute velocity
        VSYU=VYTM(j,i)+VSYR.OL(j,i);
        VSYD=VYTM(j,i)+VSYR.OL(j+1,i);
        if(VSYU>=0.0)
            VSY.OL(j,i)=VSYU;
        end
        if(VSYD<0.0)
            VSY.OL(j,i)=VSYD;
        end
        
        %OPX y-axis absolute velocity
        VSYU=VYTM(j,i)+VSYR.OPX(j,i);
        VSYD=VYTM(j,i)+VSYR.OPX(j+1,i);
        if(VSYU>=0.0)
            VSY.OPX(j,i)=VSYU;
        end
        if(VSYD<0.0)
            VSY.OPX(j,i)=VSYD;
        end
        
        %CPX y-axis absolute velocity
        VSYU=VYTM(j,i)+VSYR.CPX(j,i);
        VSYD=VYTM(j,i)+VSYR.CPX(j+1,i);
        if(VSYU>=0.0)
            VSY.CPX(j,i)=VSYU;
        end
        if(VSYD<0.0)
            VSY.CPX(j,i)=VSYD;
        end
        
        %PL y-axis absolute velocity
        VSYU=VYTM(j,i)+VSYR.PL(j,i);
        VSYD=VYTM(j,i)+VSYR.PL(j+1,i);
        if(VSYU>=0.0)
            VSY.PL(j,i)=VSYU;
        end
        if(VSYD<0.0)
            VSY.PL(j,i)=VSYD;
        end
        
        %ILM y-axis absolute velocity
        VSYU=VYTM(j,i)+VSYR.ILM(j,i);
        VSYD=VYTM(j,i)+VSYR.ILM(j+1,i);
        if(VSYU>=0.0)
            VSY.ILM(j,i)=VSYU;
        end
        if(VSYD<0.0)
            VSY.ILM(j,i)=VSYD;
        end
        
        %NOTE: VSYR.ANY are on general grid!
    end
end

for i=1:NIX+2
    %VSY.ANY(1,1:NIX+2)=0.0 --> top impermeable boudanry
    VSY.OL(1,i)=0.0;
    VSY.OPX(1,i)=0.0;
    VSY.CPX(1,i)=0.0;
    VSY.PL(1,i)=0.0;
    VSY.ILM(1,i)=0.0;
    
    %VSY.ANY(NIY+1,1:NIX+2)=0.0 --> bottom impermeable boundary
    VSY.OL(NIY+1,i)=0.0;
    VSY.OPX(NIY+1,i)=0.0;
    VSY.CPX(NIY+1,i)=0.0;
    VSY.PL(NIY+1,i)=0.0;
    VSY.ILM(NIY+1,i)=0.0;
end

% % Left + Right NO SLIP boundary
% for i=1:NIY+1
%     %VSY.ANY(1:NIY+1,1)=-VSY.ANY(1:NIY+1,2) --> left NO SLIP boundary
%     VSY.OL(i,1)=-VSY.OL(i,2);
%     VSY.OPX(i,1)=-VSY.OPX(i,2);
%     VSY.CPX(i,1)=-VSY.CPX(i,2);
%     VSY.PL(i,1)=-VSY.PL(i,2);
%     VSY.ILM(i,1)=-VSY.ILM(i,2);
%     
%     %VSY.ANY(1:NIY+1,NIX+2)=-VSY.ANY(1:NIY+1,NIX+1) --> right NO SLIP boundary
%     VSY.OL(i,NIX+2)=-VSY.OL(i,NIX+1);
%     VSY.OPX(i,NIX+2)=-VSY.OPX(i,NIX+1);
%     VSY.CPX(i,NIX+2)=-VSY.CPX(i,NIX+1);
%     VSY.PL(i,NIX+2)=-VSY.PL(i,NIX+1);
%     VSY.ILM(i,NIX+2)=-VSY.ILM(i,NIX+1);
% end

% Left + Right FREE boundary
for i=1:NIY+1
   %VSY.ANY(1:NIY+1,1)=VSY.ANY(1:NIY+1,2) --> left FREE boundary
   VSY.OL(i,1)=VSY.OL(i,2);
   VSY.OPX(i,1)=VSY.OPX(i,2);
   VSY.CPX(i,1)=VSY.CPX(i,2);
   VSY.PL(i,1)=VSY.PL(i,2);
   VSY.ILM(i,1)=VSY.ILM(i,2);

   %VSY.ANY(1:NIY+1,NIX+2)=VSY.ANY(1:NIY+1,NIX+1) --> right FREE boundary
   VSY.OL(i,NIX+2)=VSY.OL(i,NIX+1);
   VSY.OPX(i,NIX+2)=VSY.OPX(i,NIX+1);
   VSY.CPX(i,NIX+2)=VSY.CPX(i,NIX+1);
   VSY.PL(i,NIX+2)=VSY.PL(i,NIX+1);
   VSY.ILM(i,NIX+2)=VSY.ILM(i,NIX+1);
end

%Assuming solid is motionless. See also FSCR1 < 0.0 for full scale porous
%media
% VSX.OL(:,:)=0.0;
% VSX.OPX(:,:)=0.0;
% VSX.CPX(:,:)=0.0;
% VSX.PL(:,:)=0.0;
% VSX.ILM(:,:)=0.0;
% 
% VSY.OL(:,:)=0.0;
% VSY.OPX(:,:)=0.0;
% VSY.CPX(:,:)=0.0;
% VSY.PL(:,:)=0.0;
% VSY.ILM(:,:)=0.0;

end