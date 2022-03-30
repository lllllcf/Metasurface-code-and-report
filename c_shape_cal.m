% mu0 = 8.85e-12; %free space permittvity
% ep0 = 4e-7*pi; %free space permeability
% %ep0 = 1.256629e-6;
% h = 10/1000;  %高度
% w = 1/1000;  %厚度
% R = 10/1000;  %内径
% g = 0.1/1000;  %开口

mu0 = 8.85e-12; %free space permittvity
ep0 = 4e-7*pi; %free space permeability
%ep0 = 1.256629e-6;

R = 7.3/1000;  %内径
theta = 130;
g = 2*R*sind(theta/2);  %开口
w = 1.55/1000;  %厚度
h = 0.035/1000;  %高度

%增大g/R 会使得最终的reasonance 增加

Rm = R + w/2;
L = mu0*Rm*( log(8*Rm/(h+w))-0.5 );

C_gap = ep0*h*w/g + ep0*(h+w+g);
C_surf = 2*ep0*(h+w)/pi*log(4*R/g);
C = C_gap + C_surf;

f = 1/2/pi/sqrt(L*C);
disp(f/1e9);
