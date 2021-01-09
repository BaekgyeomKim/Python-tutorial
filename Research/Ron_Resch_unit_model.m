clc, clear, close all

%% Input variables
alpha_1 = 0:0.1:180; % input variable : dihedral angle between plane (13,22,23) and (13,14,24)
alpha_1_1 = 0:1:180;
alpha_2 = 120:1/30:180;

R = 25; % 20 mm % a side length of the unit triangle

%alpha_2 = 2.*atand(2.*tand(alpha_1./2)+sqrt(3.*(1+(tand(alpha_1)./2).^2))); 
% alpha_1 to alpha_2, alpha_2 represents dihedral angle between plane
% (13,14,24) and (14,15,23)

theta = 2.*asind(sqrt((1/8).*(1-cosd(alpha_1)))); 
% theta represents angle(22,13,14) in virtual plane

eta = acosd((sqrt(2.*(1-cosd(alpha_2)))-sind(theta./2))./(sqrt(3).*cosd(theta./2)));
% eta represents dihedral angle between virtual and xy plane

%% coordinates of the point 24

if sqrt(3)*tand(theta./2)>cosd(eta)
    x_24 = sqrt(3).*R.*cosd(eta).*cosd(theta./2);
else
    x_24 = R./2.*(sqrt(3).*cosd(eta).*cosd(theta./2)+3.*sind(theta./2));
end

if sqrt(3)*tand(theta./2)>cosd(eta)
    y_24 = R.*cosd(eta).*cos(theta./2);
else 
    y_24 = R./2.*(cosd(eta).*cosd(theta./2)+sqrt(3).*sind(theta./2));
end

z_24 = R.*cosd(theta./2).*sind(eta);

%% coordinates of the point 14

x_14 = R./2.*(sqrt(3).*cosd(eta).*cosd(theta./2)+sind(theta./2));
y_14 = R./2.*abs(sqrt(3).*sind(theta./2)-cosd(eta).*cosd(theta./2));
z_14 = z_24;

%% coordinates of the point 22

if sqrt(3)*tand(theta./2)>cosd(eta)
    x_22 = R.*sind(theta./2);
else
    x_22 = R./2.*(sqrt(3).*cosd(eta).*cosd(theta./2)-sind(theta./2));
end

if sqrt(3)*tand(theta./2)>cosd(eta)
    y_22 = R.*cosd(eta).*cos(theta./2);
else 
    y_22 = R./2.*(cosd(eta).*cosd(theta./2)+sqrt(3).*sind(theta./2));
end

z_22 = z_24;

%% coordinates of the point 13

x_13 = zeros(1,1801);
y_13 = zeros(1,1801);
z_13 = zeros(1,1801);

%% coordinates of the point 15

x_15 = R.*(sqrt(3).*cosd(eta).*cosd(theta./2)+sind(theta./2));
y_15 = zeros(1,181);
z_15 = zeros(1,181);

%% coordinates of the point 32

x_32 = x_15./2;
y_32 = sqrt(3)./2.*x_15;
z_32 =zeros(1,181);

%% coordinates of the point 23

x_23 = x_14;
y_23 = sqrt(3)./6.*x_15;
z_23 = -sqrt((sqrt(3)./3.*2.*R).^2-x_23.^2-y_23.^2);

%% decrease point for triangle
for i = 1:45
    x_24_1(:,i) = x_24(:,40*i); 
    y_24_1(:,i) = y_24(:,40*i); 
    x_22_1(:,i) = x_22(:,40*i); 
    y_22_1(:,i) = y_22(:,40*i);
    x_14_1(:,i) = x_14(:,40*i); 
    y_14_1(:,i) = y_14(:,40*i);
end

p_24_1 = [x_24_1;y_24_1];
p_22_1 = [x_22_1;y_22_1];
p_14_1 = [x_14_1;y_14_1];

p_24_t = transpose(p_24_1);
p_22_t = transpose(p_22_1);
p_14_t = transpose(p_14_1);

for i = 1:45
    tri(:,:,i) = [p_24_t(i,:) ; p_22_t(i,:) ; p_14_t(i,:)];
end
    

%% Estimating ellipsoid equation

a = sqrt(x_23.^2 + y_23.^2);
b = a;
c = z_14.^2.*a.^2./(x_15.*x_23 + 2.*y_14.*y_23 -(x_15.^2)./4 -(y_14.^2));

x0 = x_23;
y0 = y_23;
z0 = zeros(1,1801);

% [X,Y,Z] = ellipsoid(x0(1,10), y0(1,10), z0(1,10), a(1,10), b(1,10), c(1,10));



n = 1801;
for i = 1:n
    [A(:,:,i),B(:,:,i),C(:,:,i)] = ellipsoid(x0(1,i), y0(1,i), z0(1,i), a(1,i), b(1,i), c(1,i));
end 


%% distance from p14 to p22

d = sqrt((x_14 - x_22).^2 + (y_14 - y_22).^2);


%% maximum height point (alpha_1 = 62 deg)
% 
% p_13_max = [x_13(1,63) y_13(1,63) z_13(1,63)];
% p_14_max = [x_14(1,63) y_14(1,63) z_14(1,63)];
% p_15_max = [x_15(1,63) y_15(1,63) z_15(1,63)];
% p_22_max = [x_22(1,63) y_22(1,63) z_22(1,63)];
% p_23_max = [x_23(1,63) y_23(1,63) z_23(1,63)];
% p_24_max = [x_24(1,63) y_24(1,63) z_24(1,63)];
% p_32_max = [x_32(1,63) y_32(1,63) z_32(1,63)];
% 
% p_max = [p_13_max; p_14_max; p_15_max; p_22_max; p_23_max; p_24_max; p_32_max];
% x_max = p_max(:,1);
% y_max = p_max(:,2);
% z_max = p_max(:,3);

%% only model
Par = parula;
% In plane position of each points 
figure(1)
plot(x_14,y_14,'Color','k','Linewidth',1.5)
hold on
plot(x_22,y_22,'Color','k','Linewidth',1.5)
hold on
plot(x_24,y_24,'Color','k','Linewidth',1.5)
for i=1:45
    hold on
    patch('Vertices',tri(:,:,i),'EdgeColor',Par(i*5,:),'FaceColor','none')
end 
xlabel('Relative potision, x (mm)')
ylabel('Relative position, y (mm)')
grid on
colorbar
xlim([10 40])
ylim([0 30])
%%
% 
figure(2)
subplot(2,1,1);
plot(alpha_1,d,'Color','k','Linewidth',1.5)
ylabel('Side length of triangle (mm)')
grid on
ylim([0 30])
subplot(2,1,2);
plot(alpha_1,z_24,'Linewidth',1.5)
hold on
plot(alpha_1,z_13,'Linewidth',1.5)
hold on
xlabel('Dihedral angle, \alpha_1 (deg)')
ylabel('Out-of-plane deformation (mm)')
grid on
ylim([-5 10])
legend('Point 24 (model)','Point 13 (model)')



%% comparison between model and experiment
% 
% xyz = readmatrix('manually_tracking.xlsx');
% x13 = xyz(:,2);
% x24 = xyz(:,3);
% y13 = xyz(:,4);
% y24 = xyz(:,5);
% z13 = xyz(:,6);
% z24 = xyz(:,7);
% X24 = xyz(:,8);
% Y24 = xyz(:,9);
% Z24 = xyz(:,10);
% d_e = xyz(:,11);
% alpha_1_e = xyz(:,12);
% 
% figure(1)
% plot(alpha_1,d)
% hold on 
% plot(alpha_1_e,d_e,'--o','Color','k')
% xlabel('Dihedral angle (deg)')
% ylabel('Distance (mm)')
% grid on
% legend('Model','Experiment')
% 
% figure(2)
% plot(x_14,y_14)
% hold on
% plot(x_22,y_22)
% hold on
% plot(x_24,y_24)
% hold on
% plot(X24,Y24,'--o')
% title('XY position of single module Ron Resch')
% xlabel('x (mm)')
% ylabel('y (mm)')
% legend('P14 (model)','P22 (model)','P24 (model)','P24 (experiment)')
% grid on
% 
% figure(3)
% plot(alpha_1,z_24)
% hold on
% plot(alpha_1,z_13)
% hold on
% plot(alpha_1_e,Z24,'--o')
% xlabel('Dihedral angle (deg)')
% ylabel('z (mm)')
% grid on
% ylim([-5 10])
% legend('Z24 (model)','Z13 (model)','Z24 (experiment)')

%% 3D surface plot

figure % alpha_1 = 0
surf(A(:,:,1),B(:,:,1),C(:,:,1),'AlphaData',gradient(C(:,:,1)),'FaceAlpha','flat')
zlim([0 30])
hold on
surf(A(:,:,21),B(:,:,21),C(:,:,21),'AlphaData',gradient(C(:,:,21)),'FaceAlpha','flat')
hold on 
surf(A(:,:,41),B(:,:,41),C(:,:,41),'AlphaData',gradient(C(:,:,41)),'FaceAlpha','flat')
hold on 
surf(A(:,:,63),B(:,:,63),C(:,:,63),'AlphaData',gradient(C(:,:,63)),'FaceAlpha','flat')
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')