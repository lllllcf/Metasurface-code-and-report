%--------------------------------------------------------------------------
%   参数设置
%--------------------------------------------------------------------------
c = 3e8;                                                                    %光速
fc = 0.3e12;                                                                %频率
lambda = c/fc;                                                              %波长
R = 1000;   %m                                                              %放置信号发生源
Received = [];
t = 4;
code = [0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];
phaseChange = [0 pi];
radar_pos = code2pos(code, phaseChange, lambda, t);
%radar_pos = [11*lambda/6 9*lambda/6 7*lambda/6 5*lambda/6 3*lambda/6 lambda/6 lambda/6 3*lambda/6  5*lambda/6 7*lambda/6 9*lambda/6 11*lambda/6;
%               0   lambda/2  0 lambda/2  0 lambda/2 0   lambda/2  0 lambda/2  0 lambda/2];            %雷达振源坐标
step = 2;                                                                   %步长

for theta = 0:step:359                                                      %让目标旋转360°
    tgt_pos = [R*cosd(theta);R*sind(theta)];                                %发射机坐标
    
    %----------------------------------------------------------------------
    %   振源距离
    %----------------------------------------------------------------------
    for idx = 1:size(radar_pos,2)
        R_temp(:,idx) = norm(radar_pos(:,idx) - tgt_pos);
    end
    %----------------------------------------------------------------------
    %   接收相位
    %----------------------------------------------------------------------
    sig = sum(1.*exp(1j.*R_temp.'/lambda * 2 * pi));                             %发射信号为1到达接收时候的强度
    
    Received = [Received sig];                                              %保存输出信号
    Rx = Received;                                                     %矢量叠加
    %----------------------------------------------------------------------
    %   可视化
    %----------------------------------------------------------------------
    f1 = figure(1);
    %f1.Position = [215 361 1071 518];
    subplot(2,2,[1 3])
    plot(radar_pos(1),radar_pos(2),'bo',tgt_pos(1),tgt_pos(2),'ro');
    axis equal;
    axis([-R R -R R]);
    grid on;
    title(['角度 -> ' num2str(theta)])

    subplot(222)
    plot(0:step:theta,db(Rx));                                              %db(x)相当于20*log10(abs(x))
    axis([0 359 -inf inf]);grid on
    title('天线方向图');xlabel('角度 °');ylabel('增益 db')
    
    subplot(224)
    polarplot(deg2rad(0:step:theta),abs(Rx));
    title('极坐标天线方向图')
    drawnow()
end

function pos = code2pos(code, phase, lambda, t)
    lengthC = length(code);
    pos = zeros(2,lengthC);
    for i = 1:lengthC
        pos(1,i) = (i - (lengthC + 1)/2)*lambda/t;
        pos(2,i) = phase(code(i) + 1)/2/pi*lambda;
    end
end