clear, clc
phaseShift = [0, 60, 120, 180]./180*pi;
hold on;
axis equal;
%  rectangle('position',[0-4,0-4,8,8],'curvature',[1, 1]);
% % %rectangle('position',[1-4,0-4,8,8],'curvature',[1, 1]);
% % 
%  rectangle('position',[2-20/6,0-20/6,20/3,20/3],'curvature',[1, 1]);
% % %rectangle('position',[3-23/6,0-23/6,23/3,23/3],'curvature',[1, 1]);
% % 
%  rectangle('position',[4-8/3,0-8/3,16/3,16/3],'curvature',[1, 1]);
% % %rectangle('position',[5-11/3,0-11/3,22/3,22/3],'curvature',[1, 1]);
% % 
%  rectangle('position',[6-4/2,0-4/2,4,4],'curvature',[1, 1]);
% % %rectangle('position',[7-7/2,0-7/2,7,7],'curvature',[1, 1]);

rectangle('position',[0-4,0-4,8,8],'curvature',[1, 1]);
rectangle('position',[2-2,0-2,4,4],'curvature',[1, 1]);
rectangle('position',[4-4,0-4,8,8],'curvature',[1, 1]);
rectangle('position',[6-2,0-2,4,4],'curvature',[1, 1]);

line([0 sqrt(8)*5],[5 0]);
% line([3.5 5.5],[0 sqrt(3)*2]);
stem([3.5],[4]);
axis equal;
axis([0 7 0 10]);
