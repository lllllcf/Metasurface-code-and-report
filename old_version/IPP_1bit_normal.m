clear, clc;
D = 1; 
k = 2*3.14; 
syms theta; %仰角
syms phi; %方位角
syms f;

code = [1 0 1 0 1 0;
        0 1 0 1 0 0;
        1 0 1 0 1 0;
        0 1 0 1 0 1;
        1 0 1 0 1 0;
        0 1 0 1 0 1]; %编码 (m,n)

    
f(theta, phi) = @(theta, phi) 0; %散射函数

[r,c] = size(code);
for m = 1:r
    for n = 1:c
        f(theta, phi) = f(theta, phi) + cos(double(P(m,n,code)) + k*D*sin(theta)*((m - 1/2)*cos(phi) + (n - 1/2)*sin(phi)));
    end
end

pp;
hold on
for x = 0:0.05:3.10
    disp(x);
    for y = 0:0.05:3.10
        [X,Y,Z] = sph2cart(x,y,g(x,y,f));
        plot3([0,X],[0,Y],[0,Z]);
    end
end
%mesh(X,Y,Z);
ezmesh('x','0','z',[-20, 20, -20, 20]);
%x0 = -30:1:30;
%y0 = -30:1:30;
%mesh(x0,zeros(61),y0);

function val = g(x,y,f)
    val = double(f(x,y));
end

function res = P(m,n,code)
    if code(m,n) == 0
        res = 0;
    else
        res = pi;
    end
end


