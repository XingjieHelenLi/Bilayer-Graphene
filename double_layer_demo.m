function double_layer_demo()

%% This function is used to plot a figure to demonstrate double layer multilattice chain and interaction
%% range of a site
close all

a=[1;0];
p1=[1/2;0];
p2=[0;1];
p3=[1/2;1];

N=2;

%% Set reference configuration
ref_x0=zeros(2*N+1,2);
ref_x1=ref_x0;
ref_x2=ref_x0;
ref_x3=ref_x0;

for i=1:(2*N+1)
    
    ref_x0(i,:)=(i-N-1)*a';
    ref_x1(i,:)=ref_x0(i,:)+p1';
    ref_x2(i,:)=ref_x0(i,:)+p2';
    ref_x3(i,:)=ref_x0(i,:)+p3';
end


figure;
mrksz = 20;
lwdth = 1.5;
hold on
plot(ref_x0(:,1),ref_x0(:,2),'-bo','Markersize', 0.8*mrksz, 'Linewidth', lwdth);
plot(ref_x1(:,1),ref_x1(:,2),'-bs','Markersize', 0.8*mrksz, 'Linewidth', lwdth);
plot(ref_x2(:,1),ref_x2(:,2),'-bs','Markersize', 0.8*mrksz, 'Linewidth', lwdth);
plot(ref_x3(:,1),ref_x3(:,2),'-bo','Markersize', 0.8*mrksz, 'Linewidth', lwdth);

xi=3;
adjust1=0.15;
adjust=0.07;
lwdth2=1.5;
%% Fourteen interaction bonds
hold on
% (100)
plot([ref_x0(xi,1),ref_x0(xi+1,1)],[ref_x0(xi,2)-adjust1,ref_x0(xi+1,2)-adjust1],'k--', 'Linewidth', lwdth2);
% (001)
plot([ref_x0(xi,1),ref_x1(xi,1)], [ref_x0(xi,2)+adjust,ref_x1(xi,2)+adjust],'k--', 'Linewidth', lwdth2);
% (002)
plot([ref_x0(xi,1),ref_x2(xi,1)], [ref_x0(xi,2),ref_x2(xi,2)],'k--', 'Linewidth', lwdth2);
% (003)
plot([ref_x0(xi,1),ref_x3(xi,1)], [ref_x0(xi,2),ref_x3(xi,2)],'k--', 'Linewidth', lwdth2);
% (110)
plot([ref_x1(xi,1),ref_x0(xi+1,1)], [ref_x1(xi,2)-adjust,ref_x0(xi+1,2)-adjust],'k--', 'Linewidth', lwdth2);
% (021)
plot([ref_x1(xi,1),ref_x2(xi,1)], [ref_x1(xi,2),ref_x2(xi,2)],'k--', 'Linewidth', lwdth2);
% (023)
plot([ref_x3(xi,1),ref_x2(xi,1)], [ref_x3(xi,2)-adjust,ref_x2(xi,2)-adjust],'k--', 'Linewidth', lwdth2);
% (130)
plot([ref_x3(xi,1),ref_x0(xi+1,1)], [ref_x3(xi,2),ref_x0(xi+1,2)],'k--', 'Linewidth', lwdth2);
% (132)
plot([ref_x3(xi,1),ref_x2(xi+1,1)], [ref_x3(xi,2)+adjust,ref_x2(xi+1,2)+adjust],'k--', 'Linewidth', lwdth2);
% (031)
plot([ref_x3(xi,1),ref_x1(xi,1)], [ref_x3(xi,2),ref_x1(xi,2)],'k--', 'Linewidth', lwdth2);
% (112)
plot([ref_x2(xi+1,1),ref_x1(xi,1)], [ref_x2(xi+1,2),ref_x1(xi,2)],'k--', 'Linewidth', lwdth2);
% (122)
plot([ref_x2(xi+1,1),ref_x2(xi,1)], [ref_x2(xi+1,2)+adjust1,ref_x2(xi,2)+adjust1],'k--', 'Linewidth', lwdth2);
% (133)
plot([ref_x3(xi+1,1),ref_x3(xi,1)], [ref_x3(xi+1,2)-adjust1,ref_x3(xi,2)-adjust1],'k--', 'Linewidth', lwdth2);
% (111)
plot([ref_x1(xi+1,1),ref_x1(xi,1)], [ref_x1(xi+1,2)+adjust1,ref_x1(xi,2)+adjust1],'k--', 'Linewidth', lwdth2);


plot(ref_x0(xi,1),ref_x0(xi,2),'c.','Markersize', 2*mrksz)
% plot(ref_x1(xi,1),ref_x1(xi,2),'r.','Markersize', 1.2*mrksz)
% plot(ref_x2(xi,1),ref_x2(xi,2),'m.','Markersize', 1.2*mrksz)
% plot(ref_x3(xi,1),ref_x3(xi,2),'g.','Markersize', 1.2*mrksz)


axis([-2.1, 2.1, -0.1, 1.1])
axis equal
set(gca, 'XTick', [], 'YTick', [],'FontName', 'Times', 'FontSize', 18, 'box', 'on')

% set(gca, 'XTick', [-100], 'YTick', [-100], 'ZTick', [-100], ...
%         'FontName', 'Times', 'FontSize', 18, 'box', 'off', 'visible', 'off')






% %% Atomistic displacements
% disp_x11=zeros(2*N+1,2);
% disp_x12=disp_x11;
% disp_x21=disp_x11;
% disp_x22=disp_x11;
% 
% disp_x11=abs(rand(2*N+1,2))*0.01;
% disp_x12=abs(rand(2*N+1,2))*0.01;
% disp_x21=abs(rand(2*N+1,2))*0.01;
% disp_x22=abs(rand(2*N+1,2))*0.01;
% 
% 
% 
% 
% %% Interpolant displacements
% 
% int_disp_x11=zeros(2*N+1,2);
% int_disp_x12=disp_x11;
% int_disp_x21=disp_x11;
% int_disp_x22=disp_x11;
% 
% 
% for i=1:3:(2*N+1)-3
%    
%     int_disp_x11(i,:)=disp_x11(i,:);
%     int_disp_x12(i,:)=disp_x12(i,:);
%     int_disp_x21(i,:)=disp_x21(i,:);
%     int_disp_x22(i,:)=disp_x22(i,:);
%     
%     int_disp_x11(i+3,:)=disp_x11(i+3,:);
%     int_disp_x12(i+3,:)=disp_x12(i+3,:);
%     int_disp_x21(i+3,:)=disp_x21(i+3,:);
%     int_disp_x22(i+3,:)=disp_x22(i+3,:);
%     
%     def_grad=( disp_x11(i+3,:)- disp_x11(i,:))/3;
%     
%     p1_grad=( disp_x12(i+3,:)- disp_x12(i,:))/3;
%    
%     p2_grad=( disp_x21(i+3,:)- disp_x21(i,:))/3;
%     
%     p3_grad=( disp_x22(i+3,:)- disp_x22(i,:))/3;
%     
%     
%     for k=1:2
%     int_disp_x11(i+k,:)=disp_x11(i,:)+ k*def_grad;
%     int_disp_x12(i+k,:)=disp_x12(i,:)+k*p1_grad;
%     int_disp_x21(i+k,:)=disp_x21(i,:)+k*p2_grad;
%     int_disp_x22(i+k,:)=disp_x22(i,:)+k*p3_grad;
%     end
% end
% 
% 
% figure;hold on
% plot(ref_x11(:,1)+disp_x11(:,1), ref_x11(:,2)+disp_x11(:,2),'ro');
% plot(ref_x12(:,1)+disp_x12(:,1), ref_x12(:,2)+disp_x12(:,2),'rs');
% plot(ref_x21(:,1)+disp_x21(:,1), ref_x21(:,2)+disp_x21(:,2),'rv');
% plot(ref_x22(:,1)+disp_x22(:,1), ref_x22(:,2)+disp_x22(:,2),'rd');
% 
% plot(ref_x11(:,1)+int_disp_x11(:,1), ref_x11(:,2)+int_disp_x11(:,2),'go');
% plot(ref_x12(:,1)+int_disp_x12(:,1), ref_x12(:,2)+int_disp_x12(:,2),'gs');
% plot(ref_x21(:,1)+int_disp_x21(:,1), ref_x21(:,2)+int_disp_x21(:,2),'gv');
% plot(ref_x22(:,1)+int_disp_x22(:,1), ref_x22(:,2)+int_disp_x22(:,2),'gd');