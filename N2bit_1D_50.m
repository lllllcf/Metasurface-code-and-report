clear, clc
%全部用弧度制！
%lambda = 300;


coeT = 10; %决定间距 
phaseChange = [0 1/3*pi 2/3*pi pi];
%phaseChange = [-2*pi -4/3*pi -2/3*pi 0];
%phaseChange = linspace(0,2*pi,16);
all = 64; %几乘几
thetaAim = 30;
codeLength = 4; %编码长度, 编码长度越长，Theoretical Error越小，但是 Simulation Result 和其 Error 会变得不稳定
                    % 但编码长度小时，可能不能模拟出精度内的解。 （编码长度大会增加旁瓣）
accurate = 1; %目标精度
%amplitudeChange = [1 1 1 1];


accurate = accurate/180*pi; %遗传算法精度
maxTime = 1000;
numberCreature = 1000; %物种数
chance = 60; %交叉概率

numberPhase = length(phaseChange);
Dphi = (phaseChange(numberPhase) - phaseChange(1))/(numberPhase - 1);

%生成Amatrix
Amatrix = 1:floor(all/4);
if floor(all/4) > 6
    Amatrix = 1:6;
end

Amatrix = asin(coeT.*Dphi./2./pi./Amatrix);
usefulA = sum((abs(Amatrix) <= 1));

%均为和法线夹角
thetaAim = thetaAim/180*pi;

[sequenceH, negTheta] = generateSequence(Amatrix, thetaAim, codeLength, accurate, maxTime, numberCreature, chance, usefulA);
sequenceH = reduceUseless(sequenceH, length(Amatrix));
finalMatrix = sequenceToMatrix(sequenceH,all,numberPhase);
disp(finalMatrix);
[E2x, E2y, E2z, E3x, E3y, E3z, finalTheta] = positiveSimulation(finalMatrix,2*pi/coeT, all, phaseChange);
thetaSimulate = finalTheta;

figure(1);
surf(E3x, E3y, E3z,'EdgeColor','none');
xlabel('Ex','fontsize',12,'fontweight','b','color','r');
ylabel('Ey','fontsize',12,'fontweight','b','color','b');
zlabel('Ez','fontsize',12,'fontweight','b');
%pp;
grid off;

figure(2)
contourf(E2x, E2y, E2z,'LineStyle','none');
xlabel('sin(\theta)cos(\phi)','fontsize',12,'fontweight','b');
ylabel('sin(\theta)sin(\phi)','fontsize',12,'fontweight','b');
%pp;
axis equal

thetaAim = thetaAim/pi*180;
finalR = sprintf('Aimed Theta: %f \nTheoretical Result: %f  Error: %f \nSimulation Result: %f  Error: %f',thetaAim, negTheta, abs(thetaAim - negTheta), thetaSimulate, abs(thetaSimulate - thetaAim));
disp(finalR);


function [E2x, E2y, E2z, E3x, E3y, E3z, finalTheta] = positiveSimulation(code, KD, all,phaseChange)
    phase_change = code.*(phaseChange(length(phaseChange)) - phaseChange(1))/(length(phaseChange) - 1) + phaseChange(1);
    %amplitude_change = ; %%%%%%%%%%%%%%%%%%%%%
    %d=lmd;%the spacing distance between two neighborhood coding elements;
    stepM=500;stepN=500; %plot step
    theta=linspace(0,pi/2,stepM);
    phi=linspace(0,pi*2,stepN);
    [tt,pp]=meshgrid(theta,phi);%N*M
    E_total=zeros(stepN,stepM);
    for i=1:1:all
        for j=1:1:all
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

function res = reduceUseless(sequence, maxCode)
    for i = 1:maxCode
        if i ~= 0 && mod(i,2) == 0
            while true
                sequenceEvenP = (sequence == i);
                if sum(sequenceEvenP) > 1
                    [~, position] = sort(~sequenceEvenP);
                    position = position(1:2);
                    sequence(position(1)) = 0;
                    sequence(position(2)) = sequence(position(2))/2;
                else
                    break
                end          
            end
        end
    end
    res = sequence;
end

function [SH, negTheta] = generateSequence(Amatrix, theta, codeLength, accurate, maxTime, numberCreature, chance, usefulA)
    maxCode = length(Amatrix);
    codeSetH = unidrnd(maxCode + 1,numberCreature,codeLength) - 1;
    codeSetH(1,:) = zeros(1,codeLength);
    codeSetH = codeSetH.*round(codeSetH >= (length(Amatrix) - usefulA));
    count = 1; %次数
    while true
        %是否结束
        [cH, end0, endTheta] = checkEnd(codeSetH, theta, Amatrix, accurate, codeLength);
        if end0
            negTheta = endTheta;
            break
        end
        
        %交叉
        JCV = generateJCvector(numberCreature, chance);
        codeSetH = applyJCRule(codeSetH, JCV, numberCreature, codeLength); %考虑是否对两个set用相同的规则
        
        %变异  %暂时跳过
        
        %fitness 
        fitness = calculateFitness(codeSetH, Amatrix, codeLength, theta, numberCreature);
        
        %筛选
        [~, position] = sort(fitness);
        position = position(1:numberCreature);
        codeSetH = codeSetH(position, :);
        
        %结束
        clc;
        disp(strcat(num2str(count/maxTime*100),'%'));
        count = count + 1;
        if count == maxTime
            cH = codeSetH(1,:);
            negTheta = decode(cH ,Amatrix, codeLength)/3.14*180;
            break
        end
    end
    
    SH = cH;
end

function fitness = calculateFitness(cH, Amatrix, codeLength, theta, numberCreature)
    fitness = zeros(1, numberCreature);
    for i = 1:length(cH)
        thetaF = decode(cH(i,:), Amatrix, codeLength);
        if thetaF ~= 300
            fitness(i) = abs(thetaF - theta);
        else
            fitness(i) = 1000;
        end
    end
end

function newSet = applyJCRule(codeSet, JCV, numberCreature, codeLength)
    newSet = codeSet;
    for i = 1:numberCreature
         if JCV(i) ~= 0
            code1 = codeSet(i,:);
            code2 = codeSet(JCV(i),:);
            
            cutNumber = unidrnd(codeLength - 1);
            codePart11 = code1(1:cutNumber);
            codePart12 = code1(cutNumber + 1:end);
            codePart21 = code2(1:cutNumber);
            codePart22 = code2(cutNumber + 1:end);
            
            finalCode1 = [codePart11 codePart22];
            finalCode2 = [codePart21 codePart12];
            newSet = [newSet;finalCode1;finalCode2];
         end
    end
end


function jcv = generateJCvector(M, chance)
    jcv = zeros(1, M);
    for i = 1:M
        if unidrnd(100) <= chance
            jcv(i) = unidrnd(M);
        else
            jcv(i) = 0;
        end
    end
end

function [codeH, end0, endTheta] = checkEnd(codeSetH, theta, Amatrix, accurate, codeLength)
    codeH = zeros(codeLength);
    end0 = false;
    thetaFinal = decode(codeSetH(1,:),Amatrix,codeLength);
    endTheta = thetaFinal/3.14*180;
    if (thetaFinal ~= 300) && abs(thetaFinal - theta) < accurate 
        end0 = true;
        codeH = codeSetH(1,:);
    end
end

function angle = decode(code, Amatrix, codeLength)
    sum = 0;
    for i = 1:codeLength
        c = code(i);
        if c ~= 0
            sum = sum + sin(Amatrix(c));
        end
    end
    if (sum >= -1) && (sum <= 1)
        angle = asin(sum);
    else
        angle = 300;
    end
end

function Matrix = sequenceToMatrix(sequenceH,all,numberPhase)
    Matrix = zeros(all);
    
    arrayH = zeros(1,all);
    for i = 1:length(sequenceH)
        if sequenceH(i) ~= 0
            array = intToArray(sequenceH(i),all,numberPhase);
            for j = 1:all
                arrayH(j) = digitalAdd(arrayH(j), array(j),numberPhase);
            end
        end
    end
    
    for i = 1:all
        Matrix(i,:) = arrayH;
    end
end


function res = digitalAdd(a1,a2,numberPhase)
     res = a1 + a2;
     if a1 + a2 > numberPhase - 1
         res = res - numberPhase;
     end
end

function array = intToArray(cG, all, numberPhase)
    array = [];
    count = 0;
    if cG > 0
        dete = 0:(numberPhase - 1);
    else
        dete = (numberPhase - 1):0;
        cG = -cG;
    end
    
    while true
        for j = dete
            for i = 1:cG
                array = [array j];
                count = count + 1;
                if count == all
                    break;
                end
            end
            if count == all
                break;
            end
        end
        if count == all 
            break;
        end
    end
end