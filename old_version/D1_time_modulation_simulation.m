clear, clc
aimAngle = 70;
fc = 50e9;
c = 3e8;
lambda = c/fc;
k = 2*pi/lambda;
T = 30; %间距
inT = 300;
fp = 1e5;
Tp = 1/fp;
n = 24; %24组，共96个，最多有48个开着。
t1 = Tp./2.*(((1:n) - 1).*k.*(lambda/T + lambda/inT*3).*sind(aimAngle)./pi-1/3);
tau = Tp/3;
t1_ = t1 - Tp/4;
t2 = Tp/2 + t1;
t2_ = Tp/2 + t1_;
R = 1000;

t = 0;
interval = 1e-6;
step = 1;

Receive = zeros(1,180/step);

while true
    st = double(status(t, t1, t2, t1_, t2_, Tp, tau));
    pos = st2pos(st, T, inT, lambda, n);
    amp = reshape(st.',[1,n*4]);
    
    for theta = 0:step:180
        tgt_pos = [R*cosd(theta);R*sind(theta)];  
        R_temp = sqrt((tgt_pos(1) - pos(1,:)).^2 + (tgt_pos(2) - pos(2,:)).^2);  
        Receive(theta/step + 1) = sum(sum(amp.*exp(1j.*R_temp.'/lambda * 2 * pi)));
    end
    
    polarplot(deg2rad(0:step:180),abs(Receive));
    title(t);
    drawnow();
    t = t + interval;
end

function st = status(t, t1, t2, t1_, t2_, Tp, tau)
    t = mod(t, Tp);
    st = [((t >=t1.') & (t <= t1.' + tau)) ((t >=t2_.') & (t <= t2_.' + tau)) ((t >=t2.') & (t <= t2.' + tau)) ((t >= t1_.') & (t <= t1_.' + tau))];    
end

function pos = st2pos(st, T, inT, lambda, n)
    pos = zeros(2, n*4);
    change = (lambda./inT.*mod(4*n - 1,4) + lambda./T.*floor(n - 1))/2;
    pos(1,:) = lambda./inT.*mod((1:n*4) - 1,4) + (lambda./T + lambda./inT.*3).*floor(((1:n*4) - 1)/4) - change;
    pos(2,:) = reshape((st.*repmat([0 lambda/4 lambda/2 lambda*3/4],[n,1])).',[1,n*4]);
end
