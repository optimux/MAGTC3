function PlotTFC(TP,FSP,CLP,CSMP,RLP,VXP,VYP,PDP,Time,Step)%VSXP,VSYP,
%Plot all field variables: T, FL, CL, rho, VX, VY, P, dFS 
%suffix 'P'-->Plot
%Created 2019-10-17

%Modified for heat and species diffusion verification! 2019-11-8

global NIX
global NIY
global dxI
global dyI
global dx
global dy
global x
global y
global RL0
global CL0

%% ------------------------ Scalar Field ----------------------------------
figure(1)
subplot(2,3,1)
xx=zeros(NIX+2,1);
yy=zeros(NIY+2,1);
xx(1)=0.0;
xx(2)=0.5*dxI;
xx(3:NIX+1)=xx(2)+[1:NIX-1]*dxI;
xx(NIX+2)=xx(NIX+1)+0.5*dxI;

yy(1)=0.0;
yy(2)=0.5*dyI;
yy(3:NIY+1)=yy(2)+[1:NIY-1]*dyI;
yy(NIY+2)=yy(NIY+1)+0.5*dyI;

%T field
%TP(1,:)=TP(2,:);%For better illustration
contourf(xx,yy,TP);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
titles=['Temperature [K]      Time=',num2str(Time),' sec       Step=',num2str(Step)];
title(titles);
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;

%FST distribution
FST=FSP.OL+FSP.OPX+FSP.CPX+FSP.PL+FSP.ILM;
subplot(2,3,2)
contourf(x(1:NIX+2),y(1:NIY+2),FST*100.0);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
title('Solid Volume Fraction [Vol.%]');
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;

%MgO concentration distribution in liquid
subplot(2,3,3)
contourf(x(1:NIX+2),y(1:NIY+2),(CLP-CL0)*100.0);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
title('C_L-C_L_0 (MgO) [wt%]');
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;

%MgO old mean solid concentration
subplot(2,3,4)
contourf(x(1:NIX+2),y(1:NIY+2),CSMP*100.0);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
title('MgO Concentration (S) [wt%]');
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;

%pressure field
subplot(2,3,5)
contourf(x(1:NIX+2),y(1:NIY+2),PDP);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
title('P-P_0 [Pa]');
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;

%PL solid volume fraction
subplot(2,3,6)
contourf(x(1:NIX+2),y(1:NIY+2),FSP.PL*100.0);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
title('Plag [Vol.%]');
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;

% figure(2)
% %% ------------------------- Velocity Field -------------------------------
% %velocity field
% VXPP=zeros(NIY,NIX);
% VYPP=zeros(NIY,NIX);
% for i=1:NIX
%     for j=1:NIY
%         VXPP(j,i)=0.5*(VXP(j+1,i)+VXP(j+1,i+1));
%     end
% end
% vxmax=max(max(abs(VXPP)));
% vxmin=min(min(abs(VXPP)));
% x_scale=0.5*(vxmax+vxmin);
% 
% for i=1:NIX
%     for j=1:NIY
%         VYPP(j,i)=0.5*(VYP(j,i+1)+VYP(j+1,i+1));
%     end
% end
% vymax=max(max(abs(VYPP)));
% vymin=min(min(abs(VYPP)));
% y_scale=0.5*(vymax+vymin);
% 
% scale=sqrt(x_scale^2+y_scale^2);
% v_max=sqrt(vxmax^2+vymax^2);
% 
% refvx=zeros(1,NIX);
% refvx(3)=scale;
% refvy=zeros(1,NIX);
% VX=zeros(NIY+1,NIX);
% VX(1,:)=refvx;
% VX(2:NIY+1,:)=VXPP;
% VY=zeros(NIY+1,NIX);
% VY(1,:)=refvy;
% VY(2:NIY+1,:)=VYPP;
% 
% ya=y;
% ya(1)=0.0-0.5*dyI;
% 
% quiver(x(2:NIX+1)*1000.0,ya(1:NIY+1)*1000.0,VX,VY);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('Velocity [m/sec]');
% axis([0 x(NIX+2)*1000.0 -2.0*dyI*1000.0 y(NIY+2)*1000.0]);
% %scaletext=['Scale: ',num2str(scale),' m/s'];
% scaletext=['V_m_a_x: ',num2str(v_max,'%E'),' m/s'];
% text(2.5*dxI*1000.0,-1.2*dyI*1000.0,scaletext);

%% -------------------------- Liquid Velocity + FS -------------------------------
figure(3)
subplot(1,2,1)
%velocity field
VXPP=zeros(NIY,NIX);
VYPP=zeros(NIY,NIX);
for i=1:NIX
    for j=1:NIY
        VXPP(j,i)=0.5*(VXP(j+1,i)+VXP(j+1,i+1));
    end
end
vxmax=max(max(abs(VXPP)));
vxmin=min(min(abs(VXPP)));
x_scale=0.5*(vxmax+vxmin);

for i=1:NIX
    for j=1:NIY
        VYPP(j,i)=0.5*(VYP(j,i+1)+VYP(j+1,i+1));
    end
end
vymax=max(max(abs(VYPP)));
vymin=min(min(abs(VYPP)));
y_scale=0.5*(vymax+vymin);

scale=sqrt(x_scale^2+y_scale^2);
v_max=sqrt(vxmax^2+vymax^2);

refvx=zeros(1,NIX);
refvx(floor(NIX/2))=scale;
refvy=zeros(1,NIX);
VX=zeros(NIY+1,NIX);
VX(1,:)=refvx;
VX(2:NIY+1,:)=VXPP;
VY=zeros(NIY+1,NIX);
VY(1,:)=refvy;
VY(2:NIY+1,:)=VYPP;

ya=y;
ya(1)=0.0-0.5*dyI;
subplot(1,2,1)
contourf(x(1:NIX+2),y(1:NIY+2),TP);%FST*100.0
colorbar;
hold on
quiver(x(2:NIX+1),ya(1:NIY+1),VX,VY,'k');

axis ij
xlabel('X [m]');
ylabel('Y [m]');
title('V_L [m/sec]');
axis([0 x(NIX+2) -2.0*dyI y(NIY+2)]);
%scaletext=['Scale: ',num2str(scale),' m/s'];
scaletext=['V_m_a_x: ',num2str(v_max,'%E'),' m/s'];
text(0.5*NIX*dxI,-1.2*dyI,scaletext);

hold off
hp=colorbar;
title(hp,'T[K]');
%title(hp,'F_S [Vol.%]')

subplot(1,2,2)
contourf(x(1:NIX+2),y(1:NIY+2),RLP-RL0);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
title('\rho_L-\rho_L_0 [kg/m^3]');
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;

subplot(1,2,2)
contourf(xx,yy,TP);
axis ij
xlabel('X [m]');
ylabel('Y [m]');
titles=['Temperature [K]      Time=',num2str(Time),' sec       Step=',num2str(Step)];
title(titles);
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
%shading interp
colorbar;
hold on

[xs,ys]=meshgrid(x(2:NIX+1),y(2:NIY+1));
x0=floor(NIX/2)+1;
streamline(stream2(xs,ys,VXPP,VYPP,x(x0)*ones(NIY,1),ys(:,1)));
title('Streamlines');
axis ij
axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
hold off

drawnow
end

