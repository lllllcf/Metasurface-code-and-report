clear, clc;
t = 1;
g0 = 2;
D = 1; 
k = 2*3.14/t; 
syms theta; %仰角
syms phi; %方位角
syms f;

interval = 0.05;
endP = 3.1;



code = xlsread('code.xlsx');
f(theta, phi) = @(theta, phi) 0; %散射函数
          
code = [     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;
     0     0     0     0     1     2     3     3     0     0     1     1     3     3     3     0     1     1     2     2     0     0     0     0;];


[r,c] = size(code);
for m = 1:r
    for n = 1:c
        f(theta, phi) = f(theta, phi) + cos(double(P(m,n,code)) + k*D*sin(theta)*((m - 1/2)*cos(phi) + (n - 1/2)*sin(phi)));
    end
end

[a,b] = size(0:interval:endP);
final = zeros(b);
usex = zeros(b);
usey = zeros(b);

countx = 1;
county = 1;
pp;
hold on

for x = 0:interval:endP
    clc;
    disp(strcat(num2str(x/3.14*100),'%'));
    county = 1;
    for y = 0:interval:endP
        %[X,Y,Z] = sph2cart(x,y,g(x,y,f));
        %if (g(x,y,f) > 1)
            %plot3([0,X],[0,Y],[0,Z]);
            final(countx,county) = g(x,y,f);
            usex(countx,county) = x;
            usey(countx,county) = y;
        %end
        county = county + 1;
    end
    countx = countx + 1;
end
%ezmesh('x','0','z',[-20, 20, -20, 20]);
% 
%  figure(2);
%  pp;
%  colorbar;
%  hold on;
%  surf(flipud(final));
%  shading flat;
 
Zfinal = double((final == max(max(final)))).*final;
Xfinal = usex.*double((final == max(max(final))));
Yfinal = usey.*double((final == max(max(final))));

disp(max(max(Xfinal))/pi*180); 
disp(max(max(Yfinal))/pi*180);
% for i = 1:b
%     for j = 1:b
%         plot3([0,Xfinal(i,j)],[0,Yfinal(i,j)],[0,Zfinal(i,j)]);
%     end
% end


% figure(3);
% pp;
% colorbar;
% hold on;
% surf(flipud(A));
% shading flat;


function val = g(x,y,f)
    val = double(f(x,y));
end

function res = P(m,n,code)
    if code(m,n) == 0
        res = -pi/2;
    elseif code(m,n) == 1
        res = 0;
    elseif code(m,n) == 2
        res = pi/2;
    else
        res = pi;
    end
end