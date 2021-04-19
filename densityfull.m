function RHO=densityfull(T,P)
%For each mineral End Member, at given T and P, this function calculates
%density of each End Member;
%All End Members are listed at bottom;
%Thermoelastic data and End Members names come from "Effects of melt
%depletion on the density and seismic velocity of garnet and spinel
%lherzolite", D. L. Schutt and C. E. Lesher, 2006;
%This function is created at 2019-2-16

%T in K; P in GPa

%T=1623.5;P=0.7;

%names of End Members, 47 in total
name={'_Fo','_Fa','_OAm','_OJa','_OFs','_OMFFs','_OCaTm','_OCTTm','_OCCTm','_OMTm','_OMFTTm','_OFTm','_OCMFTm','_ODi','_OHd','_OEn','_OFes','_CAm','_CJa','_CFs','_CMFFs','_CCaTm',...
    '_CCTTm','_CCCTm','_CMTm','_CMFTTm','_CFTm','_CCMFTm','_CDi','_CHd','_CEn','_CFes','_Py','_Amd','_Sa','_Ad','_Gs','_Uv','_Sp','_Cm','_Hy','_MCm','_Mag','_Qd','_Us','_Ab','_An'};

[M,N]=size(name);
Texpan=cell(N,3);%End Member name, thermal expansivity, thermal expansivity T integration
Pexpan=cell(N,2);%End Member name, thermal expansivity P integration
RHO=cell(N,2);%End Member name, density
i=0;mem='';
ddKSdPP=0.0;%GPa^-1

%Standard T and P
P0=0.0001;%GPa
T0=298.0;%K

%% =====================End members thermoelastic data=====================

%KS: adiabatic bulk modulus,in GPa
%KT: isothermal bulk modulus,in GPa------------------------------------------
KS_Fo=128.8;
KS_Fa=138.0;

KS_OAm=106.0;
KS_OJa=124.5;
KS_OFs=124.5;
KS_OMFFs=124.5;
KS_OCaTm=207.0;
KS_OCTTm=207.0;
KS_OCCTm=207.0;
KS_OMTm=207.0;
KS_OMFTTm=207;
KS_OFTm=207.0;
KS_OCMFTm=207.0;
KS_ODi=110.5;
KS_OHd=119.0;
KS_OEn=105.6;
KS_OFes=98.8;

KS_CAm=106.0;
KS_CJa=124.5;
KS_CFs=124.5;
KS_CMFFs=124.5;
KS_CCaTm=207.0;
KS_CCTTm=207.0;
KS_CCCTm=207.0;
KS_CMTm=207.0;
KS_CMFTTm=207.0;
KS_CFTm=207.0;
KS_CCMFTm=207.0;
KS_CDi=110.5;
KS_CHd=119.0;
KS_CEn=105.6;
KS_CFes=98.8;

KS_Py=172.1;
KS_Amd=174.1;
KS_Sa=178.8;
KS_Ad=159.4;
KS_Gs=170.0;
KS_Uv=162.0;

KS_Sp=198.0;
KS_Cm=203.0;
KS_Hy=210.3;
KT_MCm=182.5;
KT_Mag=222.0;
KT_Qd=172.0;
KT_Us=222.0;

KT_Ab=57.6;
KT_An=81.6;

%dKS/dP in GPa/GPa
%dKT/dP in GPa/GPa-------------------------------------------------------------
dKSdP_Fo=4.2;
dKSdP_Fa=4.7;

dKSdP_OAm=4.8;
dKSdP_OJa=5.0;
dKSdP_OFs=5.0;
dKSdP_OMFFs=5.0;
dKSdP_OCaTm=10.2;
dKSdP_OCTTm=10.2;
dKSdP_OCCTm=10.2;
dKSdP_OMTm=10.2;
dKSdP_OMFTTm=10.2;
dKSdP_OFTm=10.2;
dKSdP_OCMFTm=10.2;
dKSdP_ODi=4.8;
dKSdP_OHd=4.0;
dKSdP_OEn=10.2;
dKSdP_OFes=9.0;

dKSdP_CAm=4.8;
dKSdP_CJa=5.0;
dKSdP_CFs=5.0;
dKSdP_CMFFs=5.0;
dKSdP_CCaTm=10.2;
dKSdP_CCTTm=10.2;
dKSdP_CCCTm=10.2;
dKSdP_CMTm=10.2;
dKSdP_CMFTTm=10.2;
dKSdP_CFTm=10.2;
dKSdP_CCMFTm=10.2;
dKSdP_CDi=4.8;
dKSdP_CHd=4.0;
dKSdP_CEn=10.2;
dKSdP_CFes=9.0;

dKSdP_Py=1.5;
dKSdP_Amd=5.85;
dKSdP_Sa=5.46;
dKSdP_Ad=3.22;
dKSdP_Gs=7.4;
dKSdP_Uv=4.7;

dKSdP_Sp=5.05;
dKSdP_Cm=5.05;
dKSdP_Hy=5.05;
dKTdP_MCm=5.8;
dKTdP_Mag=4.1;
dKTdP_Qd=4.0;
dKTdP_Us=4.1;

dKTdP_Ab=4.0;
dKTdP_An=4.0;

%dKS/dT in GPa/K
%dKT/dT in GPa/K------------------------------------------------------------
dKSdT_Fo=-0.017;
dKSdT_Fa=-0.0204;

dKSdT_OAm=-0.0205;
dKSdT_OJa=-0.0165;
dKSdT_OFs=-0.0165;
dKSdT_OMFFs=-0.0165;
dKSdT_OCaTm=-0.037;
dKSdT_OCTTm=-0.037;
dKSdT_OCCTm=-0.037;
dKSdT_OMTm=-0.037;
dKSdT_OMFTTm=-0.037;
dKSdT_OFTm=-0.037;
dKSdT_OCMFTm=-0.037;
dKSdT_ODi=-0.0205;
dKSdT_OHd=-0.0205;
dKSdT_OEn=-0.037;
dKSdT_OFes=-0.037;

dKSdT_CAm=-0.0205;
dKSdT_CJa=-0.0165;
dKSdT_CFs=-0.0165;
dKSdT_CMFFs=-0.0165;
dKSdT_CCaTm=-0.037;
dKSdT_CCTTm=-0.037;
dKSdT_CCCTm=-0.037;
dKSdT_CMTm=-0.037;
dKSdT_CMFTTm=-0.037;
dKSdT_CFTm=-0.037;
dKSdT_CCMFTm=-0.037;
dKSdT_CDi=-0.0205;
dKSdT_CHd=-0.0205;
dKSdT_CEn=-0.037;
dKSdT_CFes=-0.037;

dKSdT_Py=-0.0191;
dKSdT_Amd=-0.0204;
dKSdT_Sa=-0.0191;
dKSdT_Ad=-0.0153;
dKSdT_Gs=-0.0148;
dKSdT_Uv=-0.0148;

dKSdT_Sp=-0.015;
dKSdT_Cm=-0.015;
dKSdT_Hy=-0.015;
dKTdT_MCm=-0.02;
dKTdT_Mag=-0.02;
dKTdT_Qd=-0.02;
dKTdT_Us=-0.02;

dKTdT_Ab=-0.02;
dKTdT_An=-0.02;

%ddKSdPP in GPa^-1
%ddKTdPP in GPa^-1-------------------------------------------------------------
ddKSdP_Fo=0.0;
ddKSdP_Fa=0.0;

ddKSdP_OAm=0.0;
ddKSdP_OJa=0.0;
ddKSdP_OFs=0.0;
ddKSdP_OMFFs=0.0;
ddKSdP_OCaTm=-1.6;
ddKSdP_OCTTm=-1.6;
ddKSdP_OCCTm=-1.6;
ddKSdP_OMTm=-1.6;
ddKSdP_OMFTTm=-1.6;
ddKSdP_OFTm=-1.6;
ddKSdP_OCMFTm=-1.6;
ddKSdP_ODi=0.0;
ddKSdP_OHd=0.0;
ddKSdP_OEn=-1.6;
ddKSdP_OFes=-1.6;

ddKSdP_CAm=0.0;
ddKSdP_CJa=0.0;
ddKSdP_CFs=0.0;
ddKSdP_CMFFs=0.0;
ddKSdP_CCaTm=-1.6;
ddKSdP_CCTTm=-1.6;
ddKSdP_CCCTm=-1.6;
ddKSdP_CMTm=-1.6;
ddKSdP_CMFTTm=-1.6;
ddKSdP_CFTm=-1.6;
ddKSdP_CCMFTm=-1.6;
ddKSdP_CDi=0.0;
ddKSdP_CHd=0.0;
ddKSdP_CEn=0.0;
ddKSdP_CFes=0.0;

ddKSdP_Py=0.0;
ddKSdP_Amd=0.0;
ddKSdP_Sa=0.0;
ddKSdP_Ad=0.0;
ddKSdP_Gs=0.0;
ddKSdP_Uv=0.0;

ddKSdP_Sp=-0.65;
ddKSdP_Cm=-0.65;
ddKSdP_Hy=-0.65;
ddKTdP_MCm=0.0;%no ddKTdPP is available now, and assume as 0 to below
ddKTdP_Mag=0.0;
ddKTdP_Qd=0.0;
ddKTdP_Us=0.0;

ddKTdP_Ab=0.0;
ddKTdP_An=0.0;

%density at standard T and P, in kg/m^3------------------------------------------
r0_Fo=3230.5;
r0_Fa=4391.5;

r0_OAm=3563.0;
r0_OJa=3363.0;
r0_OFs=3363.0;
r0_OMFFs=3363.0;
r0_OCaTm=3435.0;
r0_OCTTm=3435.0;
r0_OCCTm=3435.0;
r0_OMTm=3435.0;
r0_OMFTTm=3435.0;
r0_OFTm=3435.0;
r0_OCMFTm=3435.0;
r0_ODi=3277.0;
r0_OHd=3656.0;
r0_OEn=3187.0;
r0_OFes=4004.0;

r0_CAm=3563.0;
r0_CJa=3363.0;
r0_CFs=3363.0;
r0_CMFFs=3363.0;
r0_CCaTm=3435.0;
r0_CCTTm=3435.0;
r0_CCCTm=3435.0;
r0_CMTm=3435.0;
r0_CMFTTm=3435.0;
r0_CFTm=3435.0;
r0_CCMFTm=3435.0;
r0_CDi=3277.0;
r0_CHd=3656.0;
r0_CEn=3187.0;
r0_CFes=4004.0;

r0_Py=3566.8;
r0_Amd=4316.3;
r0_Sa=4194.2;
r0_Ad=3851.3;
r0_Gs=3600.3;
r0_Uv=3847.3;

r0_Sp=3820.9;
r0_Cm=5090.0;
r0_Hy=4280.0;
r0_MCm=4414.0;
r0_Mag=5200.0;
r0_Qd=3535.0;
r0_Us=4775.0;

r0_Ab=2615.0;
r0_An=2765.0;

%alpha_0, in 10^-4---------------------------------------------------------
a0_Fo=0.285;
a0_Fa=0.2386;

a0_OAm=0.232;
a0_OJa=0.256;
a0_OFs=0.256;
a0_OMFFs=0.256;
a0_OCaTm=0.271;
a0_OCTTm=0.271;
a0_OCCTm=0.271;
a0_OMTm=0.271;
a0_OMFTTm=0.271;
a0_OFTm=0.271;
a0_OCMFTm=0.271;
a0_ODi=0.232;
a0_OHd=0.232;
a0_OEn=0.271;
a0_OFes=0.308;

a0_CAm=0.232;
a0_CJa=0.256;
a0_CFs=0.256;
a0_CMFFs=0.256;
a0_CCaTm=0.271;
a0_CCTTm=0.271;
a0_CCCTm=0.271;
a0_CMTm=0.271;
a0_CMFTTm=0.271;
a0_CFTm=0.271;
a0_CCMFTm=0.271;
a0_CDi=0.232;
a0_CHd=0.232;
a0_CEn=0.271;
a0_CFes=0.308;

a0_Py=0.2311;
a0_Amd=0.1776;
a0_Sa=0.2927;
a0_Ad=0.2103;
a0_Gs=0.1951;
a0_Uv=0.2232;

a0_Sp=0.187;
a0_Cm=0.0977;
a0_Hy=0.0513;
a0_MCm=0.143;
a0_Mag=0.635;
a0_Qd=0.2115;
a0_Us=0.635;

a0_Ab=0.0;%follow Tribaudino et al., 2010
a0_An=0.0;

%alpha_1, in 10^-8---------------------------------------------------------
a1_Fo=1.008;
a1_Fa=1.153;

a1_OAm=1.88;
a1_OJa=0.26;
a1_OFs=0.26;
a1_OMFFs=0.26;
a1_OCaTm=1.2;
a1_OCTTm=1.2;
a1_OCCTm=1.2;
a1_OMTm=1.2;
a1_OMFTTm=1.2;
a1_OFTm=1.2;
a1_OCMFTm=1.2;
a1_ODi=1.88;
a1_OHd=1.88;
a1_OEn=1.2;
a1_OFes=0.978;

a1_CAm=1.88;
a1_CJa=0.26;
a1_CFs=0.26;
a1_CMFFs=0.26;
a1_CCaTm=1.2;
a1_CCTTm=1.2;
a1_CCCTm=1.2;
a1_CMTm=1.2;
a1_CMFTTm=1.2;
a1_CFTm=1.2;
a1_CCMFTm=1.2;
a1_CDi=1.88;
a1_CHd=1.88;
a1_CEn=1.2;
a1_CFes=0.978;

a1_Py=0.5956;
a1_Amd=1.214;
a1_Sa=0.2726;
a1_Ad=0.6839;
a1_Gs=0.8089;
a1_Uv=0.5761;

a1_Sp=0.975;
a1_Cm=1.9392;
a1_Hy=1.5936;
a1_MCm=1.119;
a1_Mag=-0.7406;
a1_Qd=1.4546;
a1_Us=-0.7406;

a1_Ab=0.0;%follow Tribaudino et al., 2010
a1_An=0.0;

%alpha_2, in 1-------------------------------------------------------------
a2_Fo=-0.384;
a2_Fa=-0.0518;

a2_OAm=0.0;
a2_OJa=0.0;
a2_OFs=0.0;
a2_OMFFs=0.0;
a2_OCaTm=-0.66;
a2_OCTTm=-0.66;
a2_OCCTm=-0.66;
a2_OMTm=-0.66;
a2_OMFTTm=-0.66;
a2_OFTm=-0.66;
a2_OCMFTm=-0.66;
a2_ODi=0.0;
a2_OHd=0.0;
a2_OEn=-0.66;
a2_OFes=-0.404;

a2_CAm=0.0;
a2_CJa=0.0;
a2_CFs=0.0;
a2_CMFFs=0.0;
a2_CCaTm=-0.66;
a2_CCTTm=-0.66;
a2_CCCTm=-0.66;
a2_CMTm=-0.66;
a2_CMFTTm=-0.66;
a2_CFTm=-0.66;
a2_CCMFTm=-0.66;
a2_CDi=0.0;
a2_CHd=0.0;
a2_CEn=-0.66;
a2_CFes=-0.404;

a2_Py=-0.4538;
a2_Amd=-0.5071;
a2_Sa=-1.156;
a2_Ad=-0.2245;
a2_Gs=-0.6617;
a2_Uv=-0.2329;

a2_Sp=-0.365;
a2_Cm=0.0;
a2_Hy=0.0;
a2_MCm=-0.1063;
a2_Mag=-4.467;
a2_Qd=0.3403;
a2_Us=-4.467;

a2_Ab=0.0;%follow Tribaudino et al., 2010
a2_An=0.0;

%alpha_3, in 10^-18--------------------------------------------------------
a3_Fo=0.0;
a3_Fa=0.0;

a3_OAm=0.0;
a3_OJa=0.0;
a3_OFs=0.0;
a3_OMFFs=0.0;
a3_OCaTm=3.1;
a3_OCTTm=3.1;
a3_OCCTm=3.1;
a3_OMTm=3.1;
a3_OMFTTm=3.1;
a3_OFTm=3.1;
a3_OCMFTm=3.1;
a3_ODi=0.0;
a3_OHd=0.0;
a3_OEn=3.1;
a3_OFes=1.52;

a3_CAm=0.0;
a3_CJa=0.0;
a3_CFs=0.0;
a3_CMFFs=0.0;
a3_CCaTm=3.1;
a3_CCTTm=3.1;
a3_CCCTm=3.1;
a3_CMTm=3.1;
a3_CMFTTm=3.1;
a3_CFTm=3.1;
a3_CCMFTm=3.1;
a3_CDi=0.0;
a3_CHd=0.0;
a3_CEn=3.1;
a3_CFes=1.52;

a3_Py=0.0;
a3_Amd=0.0;
a3_Sa=0.0;
a3_Ad=0.0;
a3_Gs=0.0;
a3_Uv=0.0;

a3_Sp=0.0;
a3_Cm=0.0;
a3_Hy=0.0;
a3_MCm=0.0;%No a3 in formula, but we try to keep data as much regular as possible
a3_Mag=0.0;
a3_Qd=0.0;
a3_Us=0.0;

a3_Ab=0.0;%follow Tribaudino et al., 2010
a3_An=0.0;

%Gruneison parameter-------------------------------------------------------
ga_Fo=1.15;
ga_Fa=1.12;

ga_OAm=1.0;
ga_OJa=1.0;
ga_OFs=1.0;
ga_OMFFs=1.0;
ga_OCaTm=1.1;
ga_OCTTm=1.1;
ga_OCCTm=1.1;
ga_OMTm=1.1;
ga_OMFTTm=1.1;
ga_OFTm=1.1;
ga_OCMFTm=1.1;
ga_ODi=1.0;
ga_OHd=1.5;
ga_OEn=1.1;
ga_OFes=1.1;

ga_CAm=1.0;
ga_CJa=1.0;
ga_CFs=1.0;
ga_CMFFs=1.0;
ga_CCaTm=1.1;
ga_CCTTm=1.1;
ga_CCCTm=1.1;
ga_CMTm=1.1;
ga_CMFTTm=1.1;
ga_CFTm=1.1;
ga_CCMFTm=1.1;
ga_CDi=1.0;
ga_CHd=1.5;
ga_CEn=1.1;
ga_CFes=1.1;

ga_Py=1.29;
ga_Amd=1.29;
ga_Sa=1.29;
ga_Ad=1.38;
ga_Gs=1.38;
ga_Uv=1.38;

ga_Sp=1.1;
ga_Cm=1.1;
ga_Hy=1.1;
ga_MCm=0.0;%No gamma needed, but we try to keep data as much regular as possible
ga_Mag=0.0;
ga_Qd=0.0;
ga_Us=0.0;

ga_Ab=0.0;
ga_An=0.0;
%% ========================================================================


%% ====================End Member density at T,P===========================

%----------------------Temperature Contribution----------------------------
for i=1:N-2%for each End Member listed in Schutt et al., 2006 & Korenaga et al., 2016
    i;
    mem=['a',name{i}];%name
    Texpan{i,1}=mem;%End Member name
    %In Schutt et al., 2006, thermal expansivity is expressed as: a=a0+a1*T+a2/T^2+a3*T^4
    %while in Korenaga et al., 2016, expression is: a=a0+a1*T+a2/T. But it
    %should be: a=a0+a1*T+a2/T^2 or a=a0+a1*T+a2/T^2+a3*T^4 since a3 for
    %species in korenaga et al., 2016 is 0.
    
    p0=eval(['a0',name{i}])*10^-4;
    p1=eval(['a1',name{i}])*10^-8;
    p2=eval(['a2',name{i}]);
    p3=eval(['a3',name{i}])*10^-18;
    
    Texpan{i,2}=p0+p1*T+p2/T^2+p3*T^4;%thermal expansivity
    Texpan{i,3}=p0*(T-T0)+0.5*p1*(T^2-T0^2)-p2*(1.0/T-1.0/T0)+0.2*p3*(T^5-T0^5);%thermal expansivity T integration
    
end

% for i=N-5:N-2%for each End Member listed in 
%     mem=['a',name{i}];%name
%     Texpan{i,1}=mem;%End Member name
%     
%     p0=0.0001*eval(['a0',name{i}]);
%     p1=10^-8*eval(['a1',name{i}]);
%     p2=eval(['a2',name{i}]);
%     
%     %in Korenaga et al., 2016, thermal expansivity is expressed as:
%     %a=a0*10^-4+a1*10^-8*T+a2/T
%     Texpan{i,2}=p0+p1*T*+p2/T^2;%thermal expansivity
%     Texpan{i,3}=p0*(T-T0)+0.5*p1*(T^2-T0^2)-p2*(1.0/T-1.0/T0);%p2*log(T/T0) thermal expansivity T integration
% end

%Albite, An=0.0
b0=2.44e-5;
b1=9.0e-9;
Texpan{N-1,1}='a_Ab';
Texpan{N-1,2}=b0+2.0*b1*(T-298.0);%thermal expansivity
Texpan{N-1,3}=(b0-2.0*b1*298.0)*(T-T0)+b1*(T^2-T0^2);%thermal expansivity T integration

%Anorthite, An=100.0
b0=2.44e-5-3.1e-7+1.8e-9;
b1=9.0e-9-4.0e-11;
Texpan{N,1}='a_An';
Texpan{N,2}=b0+2.0*b1*(T-298.0);%thermal expansivity
Texpan{N,3}=(b0-2.0*b1*298.0)*(T-T0)+b1*(T^2-T0^2);%thermal expansivity T integration

%--------------------------------------------------------------------------

%-------------------------Pressure Contribution----------------------------
for i=1:N-6%for each End Member listed in Schutt et al., 2006. They use KS
    i;
    Pexpan{i,1}=['Pinteg',name{i}];
    
    %There are two regimes in ddKS/dPP
    ddKSdPP=eval(['ddKSdP',name{i}]);
    
    %Regime 1, ddKS/dPP==0.0
    if(ddKSdPP==0.0)
        A=eval(['KS',name{i}])+eval(['dKSdT',name{i}])*(T-T0);%KS(P0,T0)+dKs/dT*(TF-T0)=KS(P0,T=TF)
        B=1.0+Texpan{i,2}*T*eval(['ga',name{i}]);%1+gamma*alpha(T)*T
        dKSdP=eval(['dKSdP',name{i}]);
        Pexpan{i,2}=B*log(abs(1.0+dKSdP*(P-P0)/A))/dKSdP;
    end
    
    %Regime 2, ddKS/dPP~=0.0
    if(ddKSdPP~=0.0)
        A=eval(['KS',name{i}])+eval(['dKSdT',name{i}])*(T-T0);%KS(P0,T0)+dKSdT/dT*(TF-T0)=KS(P0,TF)
        B=1.0+Texpan{i,2}*T*eval(['ga',name{i}]);%1+gamma*alpha(T)*T
        C=eval(['dKSdP',name{i}]);%dKS/dP
        D=eval(['ddKSdP',name{i}]);%d^2KS/dP^2
        
        %integration of reciprocal quadrature
        xa=0.5*D;
        xb=C;
        xc=A-C*P0-0.5*D*P0^2;
        if(xb^2>4.0*xa*xc)
            Pexpan{i,2}=log(abs(2.0*xa*P+xb-sqrt(xb^2-4.0*xa*xc))/abs(2.0*xa*P+xb+sqrt(xb^2-4.0*xa*xc)))-log(abs(2.0*xa*P0+xb-sqrt(xb^2-4.0*xa*xc))/abs(2.0*xa*P0+xb+sqrt(xb^2-4.0*xa*xc)));
            Pexpan{i,2}=B*Pexpan{i,2}/sqrt(xb^2-4.0*xa*xc);
        else
            Pexpan{i,2}=atan((2.0*xa*P+xb)/sqrt(4.0*xa*xc-xb^2))-atan((2.0*xa*P0+xb)/sqrt(4.0*xa*xc-xb^2));
            Pexpan{i,2}=B*2.0*Pexpan{i,2}/sqrt(4.0*xa*xc-xb^2);
        end
        
    end
    
end

for i=N-5:N%for each End Member listed in Korenaga et al., 2016. They use KT
    Pexpan{i,1}=['Pinteg',name{i}];
    
    A=eval(['KT',name{i}])+eval(['dKTdT',name{i}])*(T-T0);
    B=eval(['dKTdP',name{i}]);
    Pexpan{i,2}=log(abs(1.0+B*(P-P0)/A))/B;
end

%--------------------------------------------------------------------------

%==========================================================================


%=======================End Members density at T,P=========================
for i=1:N
    RHO{i,1}=['rho',name{i}];
    rho0=eval(['r0',name{i}]);
    RHO{i,2}=rho0*exp(Pexpan{i,2}-Texpan{i,3});
end
%==========================================================================
%fprintf('End Member density returned!\n');
end

%% ==========================End-member abbreviation=======================
%Abbreviations of End Members

%Mineral                   Abbr.                    Formaula
%--------------------------olivine----------------------------------
%Forsterite                Fo                       Mg2SiO4
%Fayalite                  Fa                       Fe2SiO4

%--------------------------OPX--------------------------------------
%Acmite                    OAm                     NaFe(III)Si2O6  
%Jadeite                   OJa                     NaAlSi2O6  
%Fassaite                  OFs                     CaAl2SiO6
%(Mg,Fe)-Fassaite          OMFFs                   (Mg,Fe)Fe(III)2SiO6
%Ca-Tschermak              OCaTm                   CaAl2SiO6 
%CaTi-Tschermak            OCTTm                   CaTiAl2O6
%CrCa-Tschermak            OCCTm                   CaCr2SiO6
%Mg-Tschermak              OMTm                    MgAl2SiO6
%(Mg,Fe)Ti-Tschermak       OMFTTm                  (Mg,Fe)TiAl2O6
%Fe-Tschermak              OFTm                    Fe(II)Al2SiO6
%Cr(Mg,Fe)-Tschermak       OCMFTm                  (Mg,Fe)Cr2SiO6
%Diopside                  ODi                     CaMgSi2O6  
%Hedenbergite              OHd                     CaFeSi2O6 
%Enstatite                 OEn                     Mg2Si2O6 
%Ferrosilite               OFes                    Fe2Si2O6  

%--------------------------CPX-------------------------------------
%Acmite                    CAm                     NaFe(III)Si2O6  
%Jadeite                   CJa                     NaAlSi2O6  
%Fassaite                  CFs                     CaAl2SiO6
%(Mg,Fe)-Fassaite          CMFFs                   (Mg,Fe)Fe(III)2SiO6
%Ca-Tschermak              CCaTm                   CaAl2SiO6 
%CaTi-Tschermak            CCTTm                   CaTiAl2O6
%CrCa-Tschermak            CCCTm                   CaCr2SiO6
%Mg-Tschermak              CMTm                    MgAl2SiO6
%(Mg,Fe)Ti-Tschermak       CMFTTm                  (Mg,Fe)TiAl2O6
%Fe-Tschermak              CFTm                    Fe(II)Al2SiO6
%Cr(Mg,Fe)-Tschermak       CCMFTm                  (Mg,Fe)Cr2SiO6
%Diopside                  CDi                     CaMgSi2O6  
%Hedenbergite              CHd                     CaFeSi2O6 
%Enstatite                 CEn                     Mg2Si2O6 
%Ferrosilite               CFes                    Fe2Si2O6  

%-------------------------garnet------------------------------------
%Pyrope                    Py                      Mg3Al2[SiO4]3
%Almandine                 Amd                     Fe3Al2[SiO4]3
%Spessartite               Sa                      Mn3Al2[SiO4]3
%Andradite                 Ad                      Ca3Fe(III)2[SiO4]3
%Grossular                 Gs                      Ca3Al2[SiO4]3
%Uvarovite                 Uv                      Ca3Cr2[SiO4]3

%--------------------------spinel------------------------------------
%Spinel                    Sp                      MgAl2O4
%Chromite                  Cm                      Fe(II)Cr2O4
%Hercynite                 Hy                      Fe(II)Al2O4
%Mg-chromite               MCm                     MgCr2O4
%Magnetite                 Mag                     Fe(III)Fe(II)Fe(III)O4
%Qandilite                 Qd                      Mg2TiO4  
%Ulvospinel                Us                      TiFe(II)2O4

%--------------------------feldspar----------------------------------
%Albite                    Ab                      NaAlSi3O8
%Anorthite                 An                      CaAl2Si2O8
