clear,clc
coeT = 10; %决定间距 
phaseChange = [0 1/3*pi pi pi];

code = [2 2 2 2 2 2 2 0;
    0 0 0 0 2 2 2 2;
    2 2 2 0 0 2 0 2;
    2 2 2 2 0 0 2 0;
    0 0 2 2 2 2 0 0;
    0 0 0 0 2 2 2 2;
    2 0 0 0 0 0 2 2;
    2 2 2 2 0 0 0 0;
    0 0 2 2 2 2 0 0;
    0 0 0 2 0 2 2 2;
    2 2 0 0 0 0 2 2;
    2 2 2 0 0 0 0 0;
    0 0 2 2 2 0 0 0;
    0 0 0 2 2 2 2 2;
    2 2 0 0 0 2 0 2;
    0 2 0 2 0 2 0 0;];

all = size(code);

[E2x, E2y, E2z, E3x, E3y, E3z, finalTheta] = positiveSimulation(code,2*pi/coeT, all, phaseChange);

figure(1);
axis equal;
surf(E3x, E3y, E3z,'EdgeColor','none', 'EdgeAlpha', 0.6);
xlabel('Ex','fontsize',12,'fontweight','b','color','r');
ylabel('Ey','fontsize',12,'fontweight','b','color','b');
zlabel('Ez','fontsize',12,'fontweight','b');
pp;
axis equal;
grid off;

axis equal; 
box on;
camlight; 
lightangle(0,45);
lighting gouraud;
camproj ('perspective');
colormap jet
colorbar

figure(2)
axis equal;
contourf(E2x, E2y, E2z,'LineStyle','none');
xlabel('sin(\theta)cos(\phi)','fontsize',12,'fontweight','b');
ylabel('sin(\theta)sin(\phi)','fontsize',12,'fontweight','b');
pp;
axis equal;

disp(finalTheta);

function [E2x, E2y, E2z, E3x, E3y, E3z, finalTheta] = positiveSimulation(code, KD, all,phaseChange)
    phase_change = code.*(phaseChange(length(phaseChange)) - phaseChange(1))/(length(phaseChange) - 1) + phaseChange(1);
    amplitude_change = code; %%%%%%%%%%%%%%%%%%%%%
    %d=lmd;%the spacing distance between two neighborhood coding elements;
    stepM=500;stepN=500; %plot step
    theta=linspace(0,pi/2,stepM);
    phi=linspace(0,pi*2,stepN);
    [tt,pp]=meshgrid(theta,phi);%N*M
    E_total=zeros(stepN,stepM);
    for i=1:1:all(1)
        for j=1:1:all(2)
            E_total=E_total+ exp(-1i.*(phase_change(i,j)-KD.*((i-1/2).*cos(pp)+(j-1/2).*sin(pp)).*sin(tt)));
        end
    end
    E3x=abs(E_total).*sin(tt).*cos(pp);
    E3y=abs(E_total).*sin(tt).*sin(pp);
    E3z=abs(E_total).*cos(tt);
    
    E2x=sin(tt).*cos(pp);
    E2y=sin(tt).*sin(pp);
    E2z=abs(E_total).*cos(tt);
    
    [x0, y0] = find(abs(E3z) == max(max(abs(E3z))));
    x0 = x0(1);
    y0 = y0(1);
    finalTheta = 90 - atan(max(max(E3z))/abs(E3y(x0,y0)))/3.14*180;
end