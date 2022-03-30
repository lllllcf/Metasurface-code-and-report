lambda = 3e8/10e9;
Tin = 5;
x = -lambda/2/Tin;
y = lambda/2/Tin;
theta = 30;
phi = 1;
[res1, res2, res3] = modifyCal(x,y,theta,phi);
disp(res1);
disp(res2);
disp(res3);

% for theta = 1:90
%     for phi = 1:360
%         res = sind(theta)*cosd(phi) + sind(theta)*sind(phi) + cosd(theta);
%         if res < 0
%             disp(-1);
%         end
%     end
% end

function [res1, res2, res3] = modifyCal(x,y,theta,phi)
    a = (x-1)*sind(theta)*cosd(phi) + (y-1)*sind(theta)*sind(phi) - cosd(theta);
    b = sind(theta)*cosd(phi) + sind(theta)*sind(phi) + cosd(theta);
    res1 = abs(a) - abs(b);
    res2 = -a-b;
    res3 = b;
end