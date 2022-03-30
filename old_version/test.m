clear, clc;
% orignGraph = readmatrix('SJTU.xlsx');
orignGraph = [1 0 0 0 1 1 1 0 1 1 1;1 0 0 0 1 0 0 0 1 0 0;1 0 0 0 1 0 0 0 1 1 1;1 1 1 0 1 1 1 0 1 0 0;].*(-2) + 1;
subplot(2,2,1);
imshow(orignGraph);
title('Origin');
pp;

T = 2;
inT = 12;
arrayNum = 10; % number of antennas

fp = 1e7; % pulse frequency
fs = 3e11; %采样频率
fd = 1e7; %数据频率
fLO = 28e9; % local carrier
c = 3e8;
lambda = c/fLO;
% pulse parameters
Tp = 1/fp;
tau = Tp/3;
beta = 2*pi/lambda;

A0 = 1;

lambda = c/fLO;
z = lambda/T; % distance between two antennas

fIF=4e7;%IF Carrier Frequency fc

data=orignGraph(:);
dataLen=length(data)/fd;%data length 40*Td
t = 0:1/fs:dataLen;

dataTemp=ones((fs/fd),1)*data';%Dealing with Each Bit Td/Ts = fs/fd
dataPlot=dataTemp(:);%Makes data vector fit the time vector 扩展一下使得时间符合
t = t(1:length(dataPlot));

ModulatedData1 = A0.*cos((2*pi*(fLO)).*t)'.*dataPlot;
ModulatedData2 = A0.*cos((2*pi*(fLO)).*t)'.*dataPlot;

%% 信号生成完毕，开始传输

theta = 50; % target location
alpha = 50; % observe location

delta_t = z*sind(alpha)/c; % time delay of two signal
delta_tt = delta_t/(1/fs); % transfer the time delay to the array index 这个时间里能跑多少数据量

t11 = -Tp/6;
delta = sind(alpha)/inT/fLO;
resSig = ModulatedData1'.* pulse1(t11,tau,t11+Tp/2,Tp,t,fs,delta) - 1i.*ModulatedData2'.*pulse2(t11-Tp/4,tau,t11+Tp/4,Tp,t,fs,delta);

for i=2:arrayNum
    t1i = (((i-1)*beta*z*sind(theta))/pi-1/3)/2*Tp;
    aa = pulse1(t1i,tau,t1i+Tp/2,Tp,t,fs,delta);
    Sigi = ModulatedData1'.* pulse1(t1i,tau,t1i+Tp/2,Tp,t,fs,delta) - 1i.*ModulatedData2'.*pulse2(t1i-Tp/4,tau,t1i+Tp/4,Tp,t,fs,delta);
    resSig = resSig + circshift(Sigi, -floor(delta_tt*(i-1)));
end

RFReceived=resSig'; % received signal

%% final data
DownReceived=RFReceived.*cos((2*pi*(fLO + fp)).*t)';
RFFiltered=lowpass(DownReceived,fd,fs);

[ro, co] = size(orignGraph);
FinalData = reshape((2*RFFiltered(1:fs/fd:end)),[ro, co]);

FinalData = sign(FinalData);

subplot(2,2,2);
imshow(FinalData);
title("inT = " + int2str(inT));
pp;

%% 理想

theta = 30; % target location
alpha = 30; % observe location

delta_t = z*sind(alpha)/c; % time delay of two signal
delta_tt = delta_t/(1/fs); % transfer the time delay to the array index 这个时间里能跑多少数据量

t11 = -Tp/6;
delta = Tp/1000;
resSig = ModulatedData1'.* pulse(t11,tau,t11+Tp/2,Tp,t,fs,delta) - 1i.*ModulatedData2'.*pulse(t11-Tp/4,tau,t11+Tp/4,Tp,t,fs,delta) ;

for i=2:arrayNum
    t1i = (((i-1)*beta*z*sind(theta))/pi-1/3)/2*Tp;
    Sigi = ModulatedData1'.* pulse(t1i,tau,t1i+Tp/2,Tp,t,fs,delta) - 1i.*ModulatedData2'.*pulse(t1i-Tp/4,tau,t1i+Tp/4,Tp,t,fs,delta);
    resSig = resSig + circshift(Sigi, -floor(delta_tt*(i-1)));
end

RFReceived=resSig'; % received signal

DownReceived=RFReceived.*cos((2*pi*(fLO + fp)).*t)';
RFFiltered=lowpass(DownReceived,fd,fs);

[ro, co] = size(orignGraph);
FinalData = reshape((2*RFFiltered(1:fs/fd:end)),[ro, co]);
FinalData = sign(FinalData);

subplot(2,2,3);
imshow(FinalData);
title('Ideal');
pp;

%% Harmonic
step = 2;
p1 = [];
p_1 = [];
p5 = [];
p_7 = [];
p_3 = [];
p3 = [];
p2 = [];
p_2 = [];
ang1 = -90;
ang2 = 90;
for alpha = ang1:step:ang2
    
    delta_t = z*sind(alpha)/c; % time delay of two signal
    delta_tt = delta_t/(1/fs); % transfer the time delay to the array index 这个时间里能跑多少数据量
    t11 = -Tp/6;
    delta = Tp/1000;
    resSig = ModulatedData1'.* pulse(t11,tau,t11+Tp/2,Tp,t,fs,delta) - 1i.*ModulatedData2'.*pulse(t11-Tp/4,tau,t11+Tp/4,Tp,t,fs,delta );
    for i=2:arrayNum
        t1i = (((i-1)*beta*z*sind(theta))/pi-1/3)/2*Tp;
        Sigi = ModulatedData1'.* pulse(t1i,tau,t1i+Tp/2,Tp,t,fs,delta) - 1i.*ModulatedData2'.*pulse(t1i-Tp/4,tau,t1i+Tp/4,Tp,t,fs,delta);
        resSig = resSig + circshift(Sigi, -floor(delta_tt*(i-1)));
    end
    [f2,FFTsig]=SSBFFT(resSig,fs);
    
     p11 = find(f2>(fLO+0.7*fp),1);
     p12 = find(f2>(fLO+1.3*fp),1);
     p1 = [p1,sum(abs(FFTsig(p11:p12)))];
     
     p_11 = find(f2>(fLO-1.3*fp),1);
     p_12 = find(f2>(fLO-0.7*fp),1);
     p_1 = [p_1,sum(abs(FFTsig(p_11:p_12)))];
     
     p51 = find(f2>(fLO+4.7*fp),1);
     p52 = find(f2>(fLO+5.3*fp),1);
     p5 = [p5,sum(abs(FFTsig(p51:p52)))];
     
     p_71 = find(f2>(fLO-7.3*fp),1);
     p_72 = find(f2>(fLO-6.7*fp),1);
     p_7 = [p_7,sum(abs(FFTsig(p_71:p_72)))];
     
     p21 = find(f2>(fLO+1.7*fp),1);
     p22 = find(f2>(fLO+2.3*fp),1);
     p2 = [p2,sum(abs(FFTsig(p21:p22)))];
     
     p_21 = find(f2>(fLO-2.3*fp),1);
     p_22 = find(f2>(fLO-1.7*fp),1);
     p_2 = [p_2,sum(abs(FFTsig(p_21:p_22)))];
     
     p31 = find(f2>(fLO+2.7*fp),1);
     p32 = find(f2>(fLO+3.3*fp),1);
     p3 = [p3,sum(abs(FFTsig(p31:p32)))];
     
     p_31 = find(f2>(fLO-3.3*fp),1);
     p_32 = find(f2>(fLO-2.7*fp),1);
     p_3 = [p_3,sum(abs(FFTsig(p_31:p_32)))];
end

hold on;
subplot(2,2,4);
%plot(ang1:step:ang2, p1,ang1:step:ang2,p_1,ang1:step:ang2,p5,ang1:step:ang2,p_7);
plot(ang1:step:ang2,p_1,ang1:step:ang2,p1,ang1:step:ang2,p5,ang1:step:ang2,p_7,ang1:step:ang2,p2,ang1:step:ang2,p3,ang1:step:ang2,p_2,ang1:step:ang2,p_3);
title('Power');
pp;

%% pulse function
function res=pulse1(t1,tau,t2,Tp,totalTime,fs,delta)
    % repmat
    pp(1:floor(1.5*delta*fs)) = 0;
    pp( (floor(1.5*delta*fs) + 1 ): (floor((1.5*delta + tau)*fs)) ) = 1;
    pp( (floor((1.5*delta + tau)*fs))+1 : (floor((t2-t1 - delta/2)*fs)) ) = 0; 
    pp( (floor((t2-t1 - delta/2)*fs) + 1) : (floor((t2-t1 - delta/2 + tau)*fs))  ) = -1; 
    pp( (floor((t2-t1 - delta/2 + tau)*fs))+1 : (floor(Tp*fs)+1)) = 0;
    pp = circshift(pp, floor(t1*fs));
    pp = repmat(pp, 1, ceil(totalTime(length(totalTime))/Tp)+1);

    res = pp(1:length(totalTime));
end

function res=pulse2(t1,tau,t2,Tp,totalTime,fs,delta)
    % repmat
    pp( 1: (floor((tau)*fs)) ) = 1;
    pp( (floor((tau)*fs))+1 : (floor((t2 - t1 + 2*delta)*fs)) ) = 0;
    pp( (floor((t2 - t1 + 2*delta)*fs) + 1) : (floor((t2 - t1 + 2*delta + tau)*fs))  ) = -1;  
    pp( (floor((t2 - t1 + 2*delta + tau)*fs))+1 : (floor(Tp*fs)+1)) = 0;
    pp = circshift(pp, floor((t1 - 1.5*delta)*fs));
    pp = repmat(pp, 1, ceil(totalTime(length(totalTime))/Tp)+1);

    res = pp(1:length(totalTime));
end

%% ========================FFT Function=====================================
function [f1,FFTSSB]=SSBFFT(signal1,fs)
    L=length(signal1);
    NFFY=2^nextpow2(L);
    %To make FFT efficient, we want NFFY to be the 2's integer's power.
    FFTSSB=fft(signal1,NFFY);%to last several terms pad with zeros
    %To divide the spectrum into 2 while taking the DC term into considering 
    NumUniquePts=ceil((NFFY+1)/2);
    FFTSSB=FFTSSB(1:NumUniquePts);
    FFTSSB=abs(FFTSSB)*2/L;%*2 for SSB, /length for normalize
    %Cut the DC part and the last Part since they belongs to both sides
    FFTSSB(1)=FFTSSB(1)/2;FFTSSB(end)=FFTSSB(end)/2;
    %for SSB, divide Sampling Rate by 2
    f1=(fs/2*linspace(0,1-1/NumUniquePts,NumUniquePts))';
end

%% pulse function
function res=pulse(t1,tau,t2,Tp,totalTime,fs,delta)
    % repmat
    temp = 0:1/fs:delta;
    slope = temp/delta;
    len = length(slope);
    
    pp(1:len) = slope;
    pp(floor((delta)*fs)+1 : floor((tau)*fs)+1) = 1;
    
    
    pp(floor((tau)*fs)+1 : floor((tau)*fs)+len) = 1-slope;
    pp(floor((tau+delta)*fs)+1 : floor((t2-t1)*fs)+1) = 0;
%     % debug
%     floor((tau+delta)*fs)+1
%     floor((Tp/2)*fs)+1
%     floor((t2-t1)*fs)+1
       
    pp(floor((t2-t1)*fs)+1 : floor((t2-t1)*fs)+len) = -slope;   
    pp(floor((t2-t1+delta)*fs)+1 : floor((t2-t1+tau)*fs)+1) = -1;
       
    
    pp(floor((t2-t1+tau)*fs)+1 : floor((t2-t1+tau)*fs)+len) = slope-1;
    pp(floor((t2-t1+tau)*fs+len)+1 : (floor(Tp*fs)+1)) = 0;

%         plot(pp(1:length(pp)));
    
    pp = circshift(pp, floor(t1*fs));
    pp = repmat(pp, 1, ceil(totalTime(length(totalTime))/Tp)+1);

    res = pp(1:length(totalTime));
end