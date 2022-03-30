 clear;clc;
 f=18e+9;%frequency
 w=2*pi*f;%angular frequency
 c=3e+8;% speed of light
 lmd=c/f;%wavelength
 k=2*pi/lmd;%wave vector
 m=24;n=24;%the size of reflective plane, such as coding metasurface
codeSet = [0,1; 0,2; 0,3; 0,4; 0,5; 0,6; 0,7; 0,8; 0,9; 0,10;
1,2; 1,3; 1,4; 1,5; 1,6; 1,7; 1,8; 1,9; 1,10;
2,3;2,4;2,5;2,6;2,7;2,8;2,9;2,10;
3,4;3,5;3,6;3,7;3,8;3,9;3,10;
4,5;4,6;4,7;4,8;4,9;4,10;5,6;5,7;5,8;5,9;5,10;
6,7;6,8;6,9;6,10;7,8;7,9;7,10;8,9;8,10;9,10;];
 
 for txx = 20:20
    for i = 1:length(codeSet)
        code = codeSet(i,:);
        reflect_phi= sequenceToMatrix(code, m);%reflect phase array, user can also change the matrix according to yours’ demand.

        reflect_phi=reflect_phi.*pi;
 
        d=lmd/txx;%the spacing distance between two neighborhood coding elements;
        M=500;N=500; %plot step
        theta=linspace(0,pi/2,M);phi=linspace(0,pi*2,N);
        [tt,pp]=meshgrid(theta,phi);%N*M
        E_total=zeros(N,M);
        for ii=1:1:m
            for jj=1:1:n
                E_total=E_total+exp(-1i.*(reflect_phi(ii,jj)-k.*d.*((ii-1/2).*cos(pp)+(jj-1/2).*sin(pp)).*sin(tt)));
            end
        end
        Ex=abs(E_total).*sin(tt).*cos(pp);Ey=abs(E_total).*sin(tt).*sin(pp);Ez=abs(E_total).*cos(tt);
        Exx=sin(tt).*cos(pp);Eyy=sin(tt).*sin(pp);Ezz=abs(E_total).*cos(tt);
 
        [x0, y0] = find(abs(Ez) == max(max(abs(Ez))));
        x0 = x0(1); y0 = y0(1);
        finalTheta = 90 - atan(max(max(Ez))/abs(Ey(x0,y0)))/3.14*180;
 
        thetaStr = int2str(round(finalTheta));
        set(0,'DefaultFigureVisible', 'off');
        surf(Ex,Ey,Ez,'EdgeColor','none');
        xlabel('Ex','fontsize',12,'fontweight','b','color','r');
        ylabel('Ey','fontsize',12,'fontweight','b','color','b');
        zlabel('Ez','fontsize',12,'fontweight','b');
        grid off;
        
        folder = int2str(txx);
        smallname1 = int2str(code(1));
        smallname2 = int2str(code(2));
        name = strcat('one_bit_simulation\',folder,'\',smallname1,'-',smallname2,'--',thetaStr,'.jpg');
        %saveas(gcf, name);
    end
 end
 
function Matrix = sequenceToMatrix(sequenceH,all)
    Matrix = zeros(all);
    
    arrayH = zeros(1,all);
    for i = 1:length(sequenceH)
        if sequenceH(i) ~= 0
            array = intToArray(sequenceH(i),all);
            for j = 1:all
                arrayH(j) = digitalAdd(arrayH(j), array(j));
            end
        end
    end
    
    for i = 1:all
        Matrix(i,:) = arrayH;
    end
end

function res = digitalAdd(a1,a2)
     res = a1 + a2;
     if res == 2
         res = 0;
     end
end

function array = intToArray(cG, all)
    array = [];
    count = 0;
    if cG > 0
        dete = 0:1;
    else
        dete = 1:0;
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

