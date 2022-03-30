clear, clc
%全部用弧度制！
lambda = 300;
coeT = 1; 
coeG = 1;
D = lambda/coeT;
UN = coeG*D;
%均为和法线夹角
%sin(theta) = 　lambda/UN = lambda/(coeG*lambda/coeT) = coeT/coeG;
%   1       %% 2         3         4         5         6 %%        7         8
%   4         8         12        16        20        24
%90.0456   30.0152   19.4811   14.4849   11.5428    9.5989    8.2174    7.1844
% 24 * 24


Amatrix = [90.0456   30.0152   19.4811   14.4849   11.5428    9.5989    8.2174    7.1844]./180.*pi;
all = 24;
theta = 45/180*pi;
phi = 90/180*pi;
times = 2; %编码次数
accurate = 5/180*pi;
maxTime = 100;
M = 100; %物种数

[sequenceH, sequenceS] = generateSequence(Amatrix, theta, phi, times, accurate, maxTime, M);
finalMatrix = sequenceToMatrix(sequenceH, sequenceS,all);
disp(finalMatrix);
delete('code.xlsx');
xlswrite('code.xlsx', finalMatrix);

function [SH, SS] = generateSequence(Amatrix, theta, phi, times, accurate, maxTime, M)
    codeSetH = generateCodeSet(times, M);
    codeSetS = generateCodeSet(times, M);
    count = 1; %次数
    while true
        %是否结束
        [cH, cS, end0] = checkEnd(codeSetH, codeSetS, M, theta, phi, Amatrix, accurate, times);
        if end0
            break
        end
        
        %交叉
        JCV = generateJCvector(M);
        codeSetH = applyJCRule(codeSetH, JCV, M, times); %考虑是否对两个set用相同的规则
        codeSetS = applyJCRule(codeSetS, JCV, M, times);
        
        %变异  %暂时跳过
        
        %fitness 
        fitness = calculateFitness(codeSetH, codeSetS, Amatrix, times, theta, phi);
        
        %筛选
        [~, position] = sort(fitness);
        position = position(1:M);
        codeSetH = codeSetH(position);
        codeSetS = codeSetS(position);
        %count = count + 1;
        
        %结束
        clc;
        disp(strcat(num2str(count/maxTime*100),'%'));
        count = count + 1;
        if count == maxTime
            cH = codeSetH(1);
            cS = codeSetS(1);
            break
        end
    end
    
    SH = codeToSequence(cH, times);
    SS = codeToSequence(cS, times);
end

function fitness = calculateFitness(cH, cS, Amatrix, times, theta, phi)
    fitness = zeros(1, length(cH));
    for i = 1:length(cH)
        p1 = decode(cH(i), Amatrix, times);
        p2 = decode(cS(i), Amatrix, times);
        
        thetaF = findTheta(p1, p2, true);
        phiF = findPhi(p1, p2);
        fitness(i) = sqrt( (thetaF - theta)^2 + (phiF - phi)^2 );
    end
end

function newSet = applyJCRule(codeSet, JCV, M, times)
    newSet = codeSet;
    for i = 1:M
         if JCV(i) ~= 0
            [code1, code2] = singleJC(codeSet(i), codeSet(JCV(i)), times);
            newSet = [newSet code1 code2];
         end
    end
end

function [Rcode1, Rcode2] = singleJC(code1, code2, times)
    cutNumber = unidrnd(times - 1);
    code1 = num2str(code1);
    code2 = num2str(code2);
    codePart11 = code1(1:cutNumber);
    codePart12 = code1((cutNumber + 1):end);
    codePart21 = code2(1:cutNumber);
    codePart22 = code2((cutNumber + 1):end);
    

    Rcode1 = str2num(strcat(codePart11, codePart22));
    Rcode2 = str2num(strcat(codePart21, codePart12));
end

function jcv = generateJCvector(M)
    jcv = zeros(1, M);
    for i = 1:M
        if unidrnd(100) <= 60
            jcv(i) = unidrnd(M);
        else
            jcv(i) = 0;
        end
    end
end

function seq = codeToSequence(code, times)
    seq = [];
    code = num2str(code);
    for i = 1:times
        seq = [seq str2num(code(i))];
    end
end

function [codeH, codeS, end0] = checkEnd(codeSetH, codeSetS, M, theta, phi, Amatrix, accurate, times)
    codeH = 0;
    codeS = 0;
    end0 = false;
    for i = 1:M
        theta0 = decode(codeSetH(i), Amatrix, times);
        phi0 = decode(codeSetS(i), Amatrix, times);
        
        thetaF = findTheta(theta0, phi0, true);
        phiF = findPhi(theta0, phi0);
        if (abs(thetaF - theta) < accurate) && (abs(phiF - phi) < accurate)
            end0 = true;
            codeH = codeSetH(i);
            codeS = codeSetS(i);
            break
        end
    end
end

%function resTheta = findTheta(theta1, theta2, add)
%function resDelta = findDelta(theta1, theta2)

function angle = decode(code, Amatrix, times)
    code = num2str(code);
    angle = sameLineCalculate( Amatrix(str2num(code(1))), Amatrix(str2num(code(2))), true);
    if angle == false
        angle = 300/180*pi;
    end
    
    for i = 1:(times - 2)
        if angle == false
            angle = 300/180*pi;
            break
        end
        angle = sameLineCalculate( angle, Amatrix(str2num(code(i + 2))), true);
    end
    
    if angle == false
        angle = 300/180*pi;
    end
end

function angle = sameLineCalculate(ang1, ang2, add)
    s1 = sin(ang1);
    s2 = sin(ang2);
    if add
        if abs(s1 + s2) < 1
            angle = asin(s1 + s2);
        else
            angle = 300;
        end
    else
        if s1 - s2 < 0
            angle = 300;
        else
            angle = asin(s1 - s2);
        end
    end
end

function res = generateCodeSet(times, number)
    res = zeros(1, number);
    for i = 1:number
        for j = 1:times
            res(i) = res(i) + (unidrnd(5) + 1)*10^(j - 1);
        end
    end
end

function Matrix = sequenceToMatrix(sequenceH,sequenceS,all)
    Matrix = zeros(all);
    arrayH = sequenceToArray(sequenceH, all);
    arrayS = sequenceToArray(sequenceS, all);
    for i = 1:all
        Matrix(i,:) = arrayH;
    end
    
    for j = 1:all
        for k = 1:all
            Matrix(k,j) = digitalAdd(Matrix(k,j),arrayS(k));
        end
    end
end

function res = sequenceToArray(sequence, all)
    arrayRes = zeros(1,all);
    for i = 1:length(sequence)
        array = intToArray(sequence(i),all);
        for j = 1:all
            arrayRes(j) = digitalAdd(arrayRes(j), array(j));
        end
    end
    res = arrayRes;
end

function res = digitalAdd(a1,a2)
     res = a1 + a2;
     if a1 + a2 > 3
         res = res - 4;
     end
end

function array = intToArray(cG, all)
    array = [];
    count = 0;
    if cG > 0
        dete = 0:3;
    else
        dete = 3:0;
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


function rad = indexToRad(index,A)
    rad = A(index)/180*pi;
end


function resPhi = findPhi(theta1, theta2)
    resPhi = double(atan(sin(theta2)/sin(theta1)));
end


function resTheta = findTheta(theta1, theta2, add)
    s1 = sin(theta1) * sin(theta1);
    s2 = sin(theta2) * sin(theta2);
    if add
        res0 = sqrt(s1 + s2);
        if res0 < 1
            resTheta = double(asin(res0));
        else
            resTheta = 300;
        end
    else
        res0 = s1 - s2;
        if res0 < 0
            resTheta = 300;
        else
            res1 = sqrt(res0);
            if res1 < 1
                resTheta = double(asin(res1));
            else
                resTheta = false;
            end
        end
    end
end