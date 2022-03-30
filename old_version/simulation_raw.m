% %% m*n reflective far-field pattern
% clear;clc;
% f=18e+9;%frequency
% %w=2*pi*f;%angular frequency
% c=3e+8;% speed of light
% lmd=c/f;%wavelength
% k=2*pi/lmd;%wave vector
% m=60;n=60;%the size of reflective plane, such as coding metasurface
% %ratio_set=1; %fixed ratio between “0” and “1” coding element
% %reflect_phi=pi.*ones(m,n);%reflect phase array, user can also change the matrix according to yours’ demand.
% reflect_amp=ones(m,n);
% reflect_phi = [2 2 2 0 0 1 1 2 2 3 3 0 0 0 1];
% reflect_phi = repmat(reflect_phi, [60 4]);
% %reflect_phi=pi.*round(rand(m,n));
% % ratio_real=(m*n-sum(sum(reflect_phi./pi)))/(sum(sum(reflect_phi./pi)));
% % while ratio_real~=ratio_set 
% %     reflect_phi=pi.*round(rand(m,n));
% %     ratio_real=(m*n-sum(sum(reflect_phi./pi)))/(sum(sum(reflect_phi./pi)));
% % end
% 
% %% Generate corresponding amplitude
% code_0_amp=1;code_1_amp=1;code_2_amp=1;code_3_amp=1;%amplitude of coding elements “0” and “1”, respectively;
% for ii=1:1:m
%     for jj=1:1:n
%         if reflect_phi(ii,jj)==0
%             reflect_amp(ii,jj)=code_0_amp;
%         elseif reflect_phi(ii,jj)==1
%             reflect_amp(ii,jj)=code_1_amp;
%         elseif reflect_phi(ii,jj)==2
%             reflect_amp(ii,jj)=code_2_amp;
%         else
%             reflect_amp(ii,jj)=code_3_amp;
%         end
%     end
% end
% reflect_phi = reflect_phi.*pi/2 - pi/2;
% 
% %%
% d=lmd;%the spacing distance between two neighborhood coding elements;
% M=500;N=500; %plot step
% theta=linspace(0,pi/2,M);phi=linspace(0,pi*2,N);
% [tt,pp]=meshgrid(theta,phi);%N*M
% E_total=zeros(N,M);
% for ii=1:1:m
%     for jj=1:1:n
%         E_total=E_total+exp(1i.*(reflect_phi(ii,jj)+k.*d.*((ii-1/2).*cos(pp)+(jj-1/2).*sin(pp)).*sin(tt)));
%     end
% end
% Ex=abs(E_total).*sin(tt).*cos(pp);Ey=abs(E_total).*sin(tt).*sin(pp);Ez=abs(E_total).*cos(tt);
% Exx=sin(tt).*cos(pp);Eyy=sin(tt).*sin(pp);Ezz=abs(E_total).*cos(tt);
% 
% %%
% figure;
% surf(Ex,Ey,Ez,'EdgeColor','none');
% xlabel('Ex','fontsize',12,'fontweight','b','color','r');
% ylabel('Ey','fontsize',12,'fontweight','b','color','b');
% zlabel('Ez','fontsize',12,'fontweight','b');
% grid off;
% 
% %%
% figure(2)
% contourf(Exx,Eyy,Ezz,'LineStyle','none');
% xlabel('sin(θ)cos(φ)','fontsize',12,'fontweight','b');
% ylabel('sin(θ)sin(φ)','fontsize',12,'fontweight','b');
% axis equal

%%
 clear;clc;
 f=18e+9;%frequency
 w=2*pi*f;%angular frequency
 c=3e+8;% speed of light
 lmd=c/f;%wavelength
 k=2*pi/lmd;%wave vector
 m=24;n=24;%the size of reflective plane, such as coding metasurface
 reflect_phi= [0 1 2 3];%reflect phase array, user can also change the matrix according to yours’ demand.
 reflect_phi = repmat(reflect_phi, [24 6]);
 reflect_amp=ones(m,n);
 
 
 code_1_amp=1;code_0_amp=1;
 %amplitude of coding elements “0” and “1”, respectively;
 reflect_phi=reflect_phi.*2.*pi/3;
 
 d=lmd/2;%the spacing distance between two neighborhood coding elements;
 M=500;N=500; %plot step
 theta=linspace(0,pi/2,M);phi=linspace(0,pi*2,N);
 [tt,pp]=meshgrid(theta,phi);%N*M
 E_total=zeros(N,M);
 for ii=1:1:m
  for jj=1:1:n
  E_total=E_total+reflect_amp(ii,jj).*exp(-1i.*(reflect_phi(ii,jj)-k.*d.*((ii-1/2).*cos(pp)+(jj-1/2).*sin(pp)).*sin(tt)));
  end
 end
 Ex=abs(E_total).*sin(tt).*cos(pp);Ey=abs(E_total).*sin(tt).*sin(pp);Ez=abs(E_total).*cos(tt);
 Exx=sin(tt).*cos(pp);Eyy=sin(tt).*sin(pp);Ezz=abs(E_total).*cos(tt);
 
 figure(1);
 surf(Ex,Ey,Ez,'EdgeColor','none');
 xlabel('Ex','fontsize',12,'fontweight','b','color','r');
 ylabel('Ey','fontsize',12,'fontweight','b','color','b');
 zlabel('Ez','fontsize',12,'fontweight','b');
 grid off;
 
 figure(2)
 contourf(Exx,Eyy,Ezz,'LineStyle','none');
 xlabel('sin(θ)cos(φ)','fontsize',12,'fontweight','b');
 ylabel('sin(θ)sin(φ)','fontsize',12,'fontweight','b');
 axis equal

 [x0, y0] = find(abs(Ez) == max(max(abs(Ez))));
 x0 = x0(1); y0 = y0(1);
 finalTheta = 90 - atan(max(max(Ez))/abs(Ey(x0,y0)))/3.14*180