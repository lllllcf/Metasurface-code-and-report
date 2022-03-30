clear,clc;
%--------------------------------------------------------------------------
%   参数设置
%--------------------------------------------------------------------------
c = 3e8;                                                                    %光速
fc = 0.3e12;                                                                %频率
lambda = c/fc;                                                              %波长
R = 1000;   %m                                                              %放置信号发生源
t = 3;
code = [0 0 0 1 1 1 2 2 2 3 3 3 0 0 0 1 1 1 2 2 2 3 3 3;
0 0 0 1 1 1 2 2 2 3 3 3 0 0 0 1 1 1 2 2 2 3 3 3;
0 0 0 1 1 1 2 2 2 3 3 3 0 0 0 1 1 1 2 2 2 3 3 3;
0 0 0 1 1 1 2 2 2 3 3 3 0 0 0 1 1 1 2 2 2 3 3 3;
0 0 0 1 1 1 2 2 2 3 3 3 0 0 0 1 1 1 2 2 2 3 3 3;
0 0 0 1 1 1 2 2 2 3 3 3 0 0 0 1 1 1 2 2 2 3 3 3;
0 0 0 1 1 1 2 2 2 3 3 3 0 0 0 1 1 1 2 2 2 3 3 3;];
%code=xlsread('.xlsx'); %读取
phaseChange = [0 3/pi 2/3*pi pi];
[radar_pos_x, radar_pos_y, radar_pos_z] = code2pos(code, phaseChange, lambda, t);
%radar_pos = [11*lambda/6 9*lambda/6 7*lambda/6 5*lambda/6 3*lambda/6 lambda/6 lambda/6 3*lambda/6  5*lambda/6 7*lambda/6 9*lambda/6 11*lambda/6;
%               0   lambda/2  0 lambda/2  0 lambda/2 0   lambda/2  0 lambda/2  0 lambda/2];            %雷达振源坐标
step = 2;                                                                   %步长

[radar_pos_y,radar_pos_x] = meshgrid(radar_pos_y,radar_pos_x);
radar_pos_x = flipud(radar_pos_x);

Receive = zeros(length(0:step:360),length(0:step:90));
[tt,pp]=meshgrid(0:step:90,0:step:360);

for theta = 0:step:90
    for phi = 0:step:360
        tgt_pos = [R*sind(theta)*cosd(phi);R*sind(theta)*sind(phi);R*cosd(theta)]; 
        R_temp = sqrt((tgt_pos(1) - radar_pos_x).^2 + (tgt_pos(2) - radar_pos_y).^2 + (tgt_pos(3) - radar_pos_z).^2);
        
        Receive(phi/step + 1,theta/step + 1) = sum(sum(1.*exp(1j.*R_temp.'/lambda * 2 * pi)));
    end
end

E3x=abs(Receive).*sind(tt).*cosd(pp);
E3y=abs(Receive).*sind(tt).*sind(pp);
E3z=abs(Receive).*cosd(tt);

% surf(E3x, E3y, E3z,'EdgeColor','none');
% xlabel('Ex','fontsize',12,'fontweight','b','color','r');
% ylabel('Ey','fontsize',12,'fontweight','b','color','b');
% zlabel('Ez','fontsize',12,'fontweight','b');
% grid off;

surf(E3x,E3y,E3z,  'EdgeColor','none', 'EdgeAlpha', 0.6);
axis equal; 
box on;
camlight; 
lightangle(0,45);
lighting gouraud;
camproj ('perspective');
colormap jet
colorbar


%[r,c] = size(pp);
%xlswrite('\test.xlsx',[reshape(pp,1,r*c);reshape(tt,1,r*c);reshape(20.*log(Receive)./log(10),1,r*c)]');

function [radar_pos_x, radar_pos_y, radar_pos_z] = code2pos(code, phase, lambda, t)
    [x, y] = size(code);
    radar_pos_z = phase(code + 1)./2./pi.*lambda;
    radar_pos_y = ([1:y] - (y + 1)/2).*lambda./t;
    radar_pos_x = -([1:x] - (x + 1)/2).*lambda./t;
end