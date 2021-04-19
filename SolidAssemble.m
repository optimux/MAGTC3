function dFSSMT=SolidAssemble%(FLTM)
%This function is used to calculate solid cumulation due to movement evaluated at old time step
%Created on 2020-7-2
%Modified for upwind scheme on 2020-7-23, Tianwen-1 Mars space craft was successfully lauched in wenchang!

%======== INPUTS =======
%FSO.ANY: old time step solid volume fraction [1]
%VSXO.ANY: old time step absolute x-axis solid velocity [m/sec]
%VSYO.ANY: old time step absolute y-axis solid velocity [m/sec]

%======= OUTPUTS =======
%dFSSM.ANY: old time step solid cumulation due to solid move [1]
%dFSSMT: old time step moving solid total cumulation [1]
%dFSSMMj: old time step major element in moving solid [kg/m^3]

global NIX
global NIY
%global FSCR1
global VSXO
global VSYO
global FSO
global dtb
global dx
global dy
global dFSSM
global dFSSMMj
global RSO
global MCS
global dFSSMTe
global TCS

%% ============= Solid Cumulation =============
dFSSMT=zeros(NIY+2,NIX+2);%dFS of all solid movement
% VSXTM=VSX;
% VSYTM=VSY;
% FSX=zeros(NIY+2,NIX+1);%FS at x-axis volume faces
% for i=1:NIX+1
%     for j=1:NIY+2
%         FSX(j,i)=1.0-0.5*(FLTM(j,i)+FLTM(j,i+1));
%         if(FSX(j,i)>=FSCR1)%Solid is not moving
%             VSXTM.OL(j,i)=0.0;
%             VSXTM.OPX(j,i)=0.0;
%             VSXTM.CPX(j,i)=0.0;
%             VSXTM.PL(j,i)=0.0;
%             VSXTM.ILM(j,i)=0.0;
%         end
%     end
% end
%
% FSY=zeros(NIY+1,NIX+2);%FS at y-axis volume faces
% for i=1:NIX+2
%     for j=1:NIY+1
%         FSY(j,i)=1.0-0.5*(FLTM(j+1,i)+FLTM(j,i));
%         if(FSY(j,i)>=FSCR1)%Solid is not moving
%             VSYTM.OL(j,i)=0.0;
%             VSYTM.OPX(j,i)=0.0;
%             VSYTM.CPX(j,i)=0.0;
%             VSYTM.PL(j,i)=0.0;
%             VSYTM.ILM(j,i)=0.0;
%         end
%     end
% end

%NOTE: Indian composer Koduri Keeravaani, better known as M. M. Keeravaani, gives me a delighted music journey when I was extremely frustrated!

%FVSX.ANY=FSO.ANY*VSXO.ANY
FVSX=struct('OL',zeros(NIY+2,NIX+1),...
    'OPX',zeros(NIY+2,NIX+1),...
    'CPX',zeros(NIY+2,NIX+1),...
    'PL',zeros(NIY+2,NIX+1),...
    'ILM',zeros(NIY+2,NIX+1));

%Upwind scheme from Eq.(2-29) in Wangshanqiang
for i=1:NIX+1
    for j=1:NIY+2
        FVSX.OL(j,i)=FSO.OL(j,i)*max(VSXO.OL(j,i),0.0)+FSO.OL(j,i+1)*min(VSXO.OL(j,i),0.0);
        FVSX.OPX(j,i)=FSO.OPX(j,i)*max(VSXO.OPX(j,i),0.0)+FSO.OPX(j,i+1)*min(VSXO.OPX(j,i),0.0);
        FVSX.CPX(j,i)=FSO.CPX(j,i)*max(VSXO.CPX(j,i),0.0)+FSO.CPX(j,i+1)*min(VSXO.CPX(j,i),0.0);
        FVSX.PL(j,i)=FSO.PL(j,i)*max(VSXO.PL(j,i),0.0)+FSO.PL(j,i+1)*min(VSXO.PL(j,i),0.0);
        FVSX.ILM(j,i)=FSO.ILM(j,i)*max(VSXO.ILM(j,i),0.0)+FSO.ILM(j,i+1)*min(VSXO.ILM(j,i),0.0);
    end
end
%NOTE: if FSO(total)>=FSCR1, then VSXO.ANY==0.0 (determined in last call of NASV.m), that is, FVSX.ANY=0.0, there is no solid move in x-axis, thus no solid cumulation in x-axis direction.

%FVSY.ANY=FSO.ANY*VSYO.ANY
FVSY=struct('OL',zeros(NIY+1,NIX+2),...
    'OPX',zeros(NIY+1,NIX+2),...
    'CPX',zeros(NIY+1,NIX+2),...
    'PL',zeros(NIY+1,NIX+2),...
    'ILM',zeros(NIY+1,NIX+2));

%Upwind scheme from Eq.(2-29) in Wangshanqiang
for i=1:NIX+2
    for j=1:NIY+1
        FVSY.OL(j,i)=FSO.OL(j,i)*max(VSYO.OL(j,i),0.0)+FSO.OL(j+1,i)*min(VSYO.OL(j,i),0.0);
        FVSY.OPX(j,i)=FSO.OPX(j,i)*max(VSYO.OPX(j,i),0.0)+FSO.OPX(j+1,i)*min(VSYO.OPX(j,i),0.0);
        FVSY.CPX(j,i)=FSO.CPX(j,i)*max(VSYO.CPX(j,i),0.0)+FSO.CPX(j+1,i)*min(VSYO.CPX(j,i),0.0);
        FVSY.PL(j,i)=FSO.PL(j,i)*max(VSYO.PL(j,i),0.0)+FSO.PL(j+1,i)*min(VSYO.PL(j,i),0.0);
        FVSY.ILM(j,i)=FSO.ILM(j,i)*max(VSYO.ILM(j,i),0.0)+FSO.ILM(j+1,i)*min(VSYO.ILM(j,i),0.0);
    end
end
%NOTE: if FSO(total)>=FSCR1, then VSYO.ANY==0.0 (determined in last call of NASV.m), that is, FVSY.ANY=0.0, there is no solid move in y-axis, thus no solid cumulation in y-axis direction.

%dFSSM is a global struct defined in MAGTC3.m
%dFSSM=struct('OL',zeros(NIY+2,NIX+2),...
%     'OPX',zeros(NIY+2,NIX+2),...
%     'CPX',zeros(NIY+2,NIX+2),...
%     'PL',zeros(NIY+2,NIX+2),...
%     'ILM',zeros(NIY+2,NIX+2));%dFS of each solid movement

for i=1:NIX
    for j=1:NIY
        dFSSM.OL(j+1,i+1)=dtb*(FVSX.OL(j+1,i+1)-FVSX.OL(j+1,i))/dx(i)+dtb*(FVSY.OL(j+1,i+1)-FVSY.OL(j,i+1))/dy(j);
        dFSSM.OPX(j+1,i+1)=dtb*(FVSX.OPX(j+1,i+1)-FVSX.OPX(j+1,i))/dx(i)+dtb*(FVSY.OPX(j+1,i+1)-FVSY.OPX(j,i+1))/dy(j);
        dFSSM.CPX(j+1,i+1)=dtb*(FVSX.CPX(j+1,i+1)-FVSX.CPX(j+1,i))/dx(i)+dtb*(FVSY.CPX(j+1,i+1)-FVSY.CPX(j,i+1))/dy(j);
        dFSSM.PL(j+1,i+1)=dtb*(FVSX.PL(j+1,i+1)-FVSX.PL(j+1,i))/dx(i)+dtb*(FVSY.PL(j+1,i+1)-FVSY.PL(j,i+1))/dy(j);
        dFSSM.ILM(j+1,i+1)=dtb*(FVSX.ILM(j+1,i+1)-FVSX.ILM(j+1,i))/dx(i)+dtb*(FVSY.ILM(j+1,i+1)-FVSY.ILM(j,i+1))/dy(j);
    end
end

%To make sure that FSO.ANY(1,2:NIX+1)=FSO.ANY(2,2:NIX+1), FSO.ANY(NIY+2,2:NIX+1)=FSO.ANY(NIY+1,2:NIX+1), FSO.ANY(1:NIY+2,1)=FSO.ANY(1:NIY+2,2),
%FSO.ANY(1:NIY+2,NIX+2)=FSO.ANY(1:NIY+2,NIX+1)
dFSSM.OL(1,2:NIX+1)=dFSSM.OL(2,2:NIX+1);
dFSSM.OL(NIY+2,2:NIX+1)=dFSSM.OL(NIY+1,2:NIX+1);
dFSSM.OL(1:NIY+2,1)=dFSSM.OL(1:NIY+2,2);
dFSSM.OL(1:NIY+2,NIX+2)=dFSSM.OL(1:NIY+2,NIX+1);

dFSSM.OPX(1,2:NIX+1)=dFSSM.OPX(2,2:NIX+1);
dFSSM.OPX(NIY+2,2:NIX+1)=dFSSM.OPX(NIY+1,2:NIX+1);
dFSSM.OPX(1:NIY+2,1)=dFSSM.OPX(1:NIY+2,2);
dFSSM.OPX(1:NIY+2,NIX+2)=dFSSM.OPX(1:NIY+2,NIX+1);

dFSSM.CPX(1,2:NIX+1)=dFSSM.CPX(2,2:NIX+1);
dFSSM.CPX(NIY+2,2:NIX+1)=dFSSM.CPX(NIY+1,2:NIX+1);
dFSSM.CPX(1:NIY+2,1)=dFSSM.CPX(1:NIY+2,2);
dFSSM.CPX(1:NIY+2,NIX+2)=dFSSM.CPX(1:NIY+2,NIX+1);

dFSSM.PL(1,2:NIX+1)=dFSSM.PL(2,2:NIX+1);
dFSSM.PL(NIY+2,2:NIX+1)=dFSSM.PL(NIY+1,2:NIX+1);
dFSSM.PL(1:NIY+2,1)=dFSSM.PL(1:NIY+2,2);
dFSSM.PL(1:NIY+2,NIX+2)=dFSSM.PL(1:NIY+2,NIX+1);

dFSSM.ILM(1,2:NIX+1)=dFSSM.ILM(2,2:NIX+1);
dFSSM.ILM(NIY+2,2:NIX+1)=dFSSM.ILM(NIY+1,2:NIX+1);
dFSSM.ILM(1:NIY+2,1)=dFSSM.ILM(1:NIY+2,2);
dFSSM.ILM(1:NIY+2,NIX+2)=dFSSM.ILM(1:NIY+2,NIX+1);

%Total dFSSM
dFSSMT=dFSSM.OL+dFSSM.OPX+dFSSM.CPX+dFSSM.PL+dFSSM.ILM;

%% ========== Major Element in Moving Solid =========
%NOTE: RSO comes from updates RS in MAGTC3.m.
%NOTE: MCS comes from updates in MAGTFC.m.
for i=1:NIX+2
    for j=1:NIY+2
        dFSSMMj.SiO2(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.SiO2.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.SiO2.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.SiO2.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.SiO2.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.SiO2.ILM(j,i);
        dFSSMMj.TiO2(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.TiO2.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.TiO2.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.TiO2.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.TiO2.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.TiO2.ILM(j,i);
        dFSSMMj.Al2O3(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.Al2O3.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.Al2O3.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.Al2O3.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.Al2O3.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.Al2O3.ILM(j,i);
        dFSSMMj.FeO(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.FeO.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.FeO.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.FeO.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.FeO.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.FeO.ILM(j,i);
        dFSSMMj.Fe2O3(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.Fe2O3.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.Fe2O3.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.Fe2O3.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.Fe2O3.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.Fe2O3.ILM(j,i);
        dFSSMMj.MnO(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.MnO.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.MnO.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.MnO.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.MnO.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.MnO.ILM(j,i);
        dFSSMMj.MgO(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.MgO.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.MgO.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.MgO.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.MgO.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.MgO.ILM(j,i);
        dFSSMMj.CaO(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.CaO.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.CaO.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.CaO.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.CaO.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.CaO.ILM(j,i);
        dFSSMMj.Na2O(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.Na2O.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.Na2O.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.Na2O.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.Na2O.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.Na2O.ILM(j,i);
        dFSSMMj.K2O(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.K2O.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.K2O.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.K2O.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.K2O.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.K2O.ILM(j,i);
        dFSSMMj.P2O5(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.P2O5.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.P2O5.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.P2O5.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.P2O5.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.P2O5.ILM(j,i);
        dFSSMMj.H2O(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*MCS.H2O.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*MCS.H2O.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*MCS.H2O.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*MCS.H2O.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*MCS.H2O.ILM(j,i);
        
        dFSSMTe.Sm(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*TCS.Sm.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*TCS.Sm.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*TCS.Sm.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*TCS.Sm.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*TCS.Sm.ILM(j,i);
        dFSSMTe.Nd(j,i)=dFSSM.OL(j,i)*RSO.OL(j,i)*TCS.Nd.OL(j,i)+dFSSM.OPX(j,i)*RSO.OPX(j,i)*TCS.Nd.OPX(j,i)+dFSSM.CPX(j,i)*RSO.CPX(j,i)*TCS.Nd.CPX(j,i)+dFSSM.PL(j,i)*RSO.PL(j,i)*TCS.Nd.PL(j,i)+dFSSM.ILM(j,i)*RSO.ILM(j,i)*TCS.Nd.ILM(j,i);
        
    end
end

end

