function [] = pp()
    fontSize = 16; %统一字号
    axis = 1;%坐标轴粗细
    delta = 4; %坐标轴数字和其余大小的差异
    set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);%线宽
    set(gca,'linewidth',axis); %坐标轴粗细
    set(gca,'FontSize',fontSize - delta);  %改变图中坐标的字号
    set(get(gca,'XLabel'),'FontSize',fontSize); %xlabel字号
    set(get(gca,'YLabel'),'FontSize',fontSize); %ylabel轴字号
    set(get(gca,'title'),'FontSize',fontSize); %title轴字号
    set(get(gca,'legend'),'FontSize',fontSize); %legend轴字号
end
