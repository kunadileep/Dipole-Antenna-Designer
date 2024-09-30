% MATLAB
clear all;
clc;

l = input("Enter the length (in m): ");
w = input("Enter the width (in m): ");

fc = input("Enter the frequency (in Hz): ");
Im = input("Enter the mximum current phasor: ");

c = 299792458;
n = 377;
y = 0.57721566;

lamda = c/fc;

%% Impedance
a = w/4;

f = (9*fc)/10:fc/100: 11*fc/10;
kL = (2*pi*l/c).*f;

p = (n./(2.*pi.*(sin(kL./2)).^2));

Rp1 = y+log(kL)-cosint(kL);
Rp2 = ((0.5).*sin(kL).*(sinint(2.*kL)-(2.*sinint(kL))));
Rp3 = ((0.5).*cos(kL).*(cosint(2.*(kL))-(2.*cosint(kL))+y+log((0.5).*kL)));
R = p.*(Rp1+Rp2+Rp3);

Xp1 = sinint(kL);
Xp2 = (0.5).*cos(kL).*(sinint(2.*kL)-(2.*sinint(kL))).*(-1);
Xp3 = (0.5).*sin(kL).*(cosint(2.*(kL))-(2.*cosint(kL))+ cosint(2*((a/l)^2)).*kL);
X = p.*(Xp1+Xp2+Xp3);

figure(1);
plot(f, R, 'b', 'DisplayName', 'Resistance (R)');
title('Impedance');
xlabel('Frequency (Hz)');
ylabel('Impedance (ohm)');
hold on;
plot(f, X, 'r', 'DisplayName', 'Reactance (X)');
hold 'off';
legend;

%% Current Distribution

z = -l/2:l/1000:l/2;
Iz = zeros(length(z));
b = (2*pi)/lamda;

for i = 1:length(z)
    if(z(i)<0)
          Iz(i) = Im*sin(b*((l/2)+z(i)));
    else
        Iz(i) = Im*sin(b*((l/2)-z(i)));
    end
end
figure(3);

plot(z,Iz,"LineWidth",2)
title("Current Distribution PLot (Iz vs z)")
xlabel('z(m)');
ylabel('Iz(A)');
%% Radiation patterns 3D and 2D

%Defining variables in spherical coordinates
theta=0:0.05:2*pi;%theta vector
phi=0:0.05:2*pi;%phi vector

l_lamda1=l/lamda;% length of antenna in terms of wavelengths 

% evaluating radiation intensity(U)
U1=((cos(l_lamda1*pi*cos(theta-(pi/2))) - cos(l_lamda1*pi) )./ sin(theta-(pi/2)) ).^2 ;

%converting to dB scale
U1_1=10*log10(U1); % add a very small number for l = even multiples of lamda.

%normalizing in order to make U vector positive
min1=min(U1_1);
U=U1_1-min1;

figure(4);
polarplot(theta, U)
title("Elevation Plot (in dBi)");
figure(5);
polarplot(phi,zeros(length(phi))+U(1,1))
title("Azimuthal Plot (in dBi)");


U(1,1)=0;
for n=1:length(phi)
    theta(n,:)=theta(1,:);
end

phi=phi'; % transpose;
for m=1:length(phi)
    phi(:,m)=phi(:,1);
end

for k=1:length(U)
    U(k,:)=U(1,:);
end

[k,i,m]=sph2cart(phi,theta,U);

%plotting
figure(6);
grid off;
C=sqrt(k.^2+i.^2+m.^2);
surf(k, i, m,C,EdgeColor="none")
colormap('jet');
o = colorbar("east");
o.Label.String = "dBi";
o.Label.HorizontalAlignment = "right";
title('Radiation Power Pattern for Dipole Antenna')
xlabel('x');
ylabel('y');
zlabel('z');
shading interp;
            
%% physical image
figure(7);
L = l; 
W = w; 
k = [-W/2, -W/2; W/2, W/2];
z = [-L/2, L/2; -L/2, L/2];
x = zeros(size(k));
mesh(x, k, z, 'EdgeColor', 'k', 'FaceColor', [0.9290 0.6940 0.1250]);
title('Dipole Strip Antenna Structure');
ylabel('y(m)');
zlabel('z(m)');
axis('equal');
axis('vis3d');
xlim([-0.001, 0.001]);
