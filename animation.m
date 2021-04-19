%This is used to generate animation
%2021-3-28

W=20.0;%width [m]
H=20.0;%height [m]

dir='20214142026';

Tfile=[pwd,'\',dir,'\TRD.mat'];
tfile=[pwd,'\',dir,'\Time.mat'];
VXfile=[pwd,'\',dir,'\VXRD.mat'];
VYfile=[pwd,'\',dir,'\VYRD.mat'];
load(Tfile);
load(tfile);
load(VXfile);
load(VYfile);

[NY,NX,frame]=size(TRD);
NIX=NX-2;%number of unit, x-axis
NIY=NY-2;%number of unit, y-axis
dxI=W/NIX;%x-axis interval [m]
dyI=H/NIY;%y-axis interval [m]

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

VXPP=zeros(NIY,NIX);
VYPP=zeros(NIY,NIX);

ya=y;
ya(1)=0.0-0.5*dyI;
refvx=zeros(1,NIX);
refvy=zeros(1,NIX);
VX=zeros(NIY+1,NIX);
VY=zeros(NIY+1,NIX);

anifile=[pwd,'\',dir,'\demo13TV.avi'];
aviobj=VideoWriter(anifile,'Uncompressed AVI');
open(aviobj);

for k=1:frame
    
    for i=1:NIX
        for j=1:NIY
            VXPP(j,i)=0.5*(VXRD(j+1,i,k)+VXRD(j+1,i+1,k));
        end
    end
    vxmax=max(max(abs(VXPP)));
    vxmin=min(min(abs(VXPP)));
    x_scale=0.5*(vxmax+vxmin);
    
    for i=1:NIX
        for j=1:NIY
            VYPP(j,i)=0.5*(VYRD(j,i+1,k)+VYRD(j+1,i+1,k));
        end
    end
    vymax=max(max(abs(VYPP)));
    vymin=min(min(abs(VYPP)));
    y_scale=0.5*(vymax+vymin);
    
    scale=sqrt(x_scale^2+y_scale^2);
    v_max=sqrt(vxmax^2+vymax^2);
    
    refvx(floor(NIX/2))=scale;
    VX(1,:)=refvx;
    VX(2:NIY+1,:)=VXPP;
    VY(1,:)=refvy;
    VY(2:NIY+1,:)=VYPP;
    
    contourf(x,y,TRD(:,:,k));
    colorbar;
    hold on
    quiver(x(2:NIX+1),ya(1:NIY+1),VX,VY,'k');
    
    axis ij
    xlabel('X [m]');
    ylabel('Y [m]');
    title(['V_L [m/sec]      Time=',num2str(t(k)),'sec']);
    axis([0 x(NIX+2) -2.0*dyI y(NIY+2)]);
    %scaletext=['Scale: ',num2str(scale),' m/s'];
    scaletext=['V_m_a_x: ',num2str(v_max,'%E'),' m/s'];
    text(0.5*NIX*dxI,-1.2*dyI,scaletext);
    
    hold off
    hp=colorbar;
    title(hp,'T[K]');
    Current_Frame=getframe(gcf);
    writeVideo(aviobj,Current_Frame);
    
end
close(aviobj);

RLfile=[pwd,'\',dir,'\RLRD.mat'];
load(RLfile);

anifile=[pwd,'\',dir,'\demo13TRL.avi'];
aviobj=VideoWriter(anifile,'Uncompressed AVI');
open(aviobj);

for i=1:frame

    contourf(x,y,RLRD(:,:,i));
    axis ij
    xlabel('X [m]');
    ylabel('Y [m]');
    titles=['\rho_l [kg.m^-^3]      Time=',num2str(t(i)),' sec'];
    title(titles);
    axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
    %shading interp
    colorbar;

    Current_Frame=getframe(gcf);
    writeVideo(aviobj,Current_Frame);

end
close(aviobj);

MCLfile=[pwd,'\',dir,'\MCLRD.mat'];
load(MCLfile);

anifile=[pwd,'\',dir,'\demo13TCL.avi'];
aviobj=VideoWriter(anifile,'Uncompressed AVI');
open(aviobj);

for i=1:frame

    contourf(x,y,MCLRD.Al2O3(:,:,i)*100.0-7.33);
    axis ij
    xlabel('X [m]');
    ylabel('Y [m]');
    titles=['MgO [wt%]      Time=',num2str(t(i)),' sec'];
    title(titles);
    axis([0.0 x(NIX+2) 0.0 y(NIY+2)]);
    %shading interp
    colorbar;

    Current_Frame=getframe(gcf);
    writeVideo(aviobj,Current_Frame);

end
close(aviobj);