% Range of frequencies to plot
% In the real index of refraction and coefficient of extinction plots
plot_range = [10^15, 10^16];

% Dielectric Properties
fs = [0.5, 0.5]; % The fraction of electrons that have the corresponding resonant frequency
ws = [3.5*10^15, 7*10^(15)]; % The electron resonant frequencies
gs = [10^10, 10^10]; % The electron damping factors

Nb = 10^27; % The electron number density

% The index of refraction of the material
% from which the incident light is coming
n_2 = 1;

% Metalic Properties
wp = 0;
%wp = 10^15 * 2 * pi;
tau = 10^(-13) * 6;

% Constants
c = 3*10^8;
e = 1.6*10^(-19);
m = 9*10^(-31);
e0 = 8.9 * 10^(-12);

% Amplitude of electron current oscillations
A = Nb*e^2/(m*e0);

% The complex index of refraction
syms w;
N_squared = ( 1 + A*sum(fs./(ws.^2-w^2-1i*w*gs)) - wp^2/(w^2+1i*w/tau) );
N = sqrt(N_squared);

% Prepare the plot
tiledlayout(2,2)

%% Index of refraction
nexttile
fplot(real(N),plot_range)
ylabel('n')
xlabel('\omega')
%ylim([0.85,1.5])
%yticks(1)
xticks(ws)
xticklabels(arrayfun(@(i) append("\omega_{",string(i),"}"), 1:length(ws)))
%xticks([wp/2 wp 1.5*wp])
%xticklabels({'0.5\omega_p' '\omega_p' '1.5\omega_p'})


%% Coefficient of extinction
nexttile
fplot(imag(N),plot_range)
ylabel('\kappa')
xlabel('\omega')
%ylim([0.85,1.5])
%yticks(1)
xticks(ws)
xticklabels(arrayfun(@(i) append("\omega_{",string(i),"}"), 1:length(ws)))
%xticks([wp/2, wp, 1.5*wp])
%xticklabels({'0.5\omega_p', '\omega_p', '1.5\omega_p'})

%% Reflectance
nexttile

% The frequency of light for which to plot the reflectance

%w = wp;
w = ws(1);

%title('\omega = \omega_p')
title('\omega = \omega_1')

theta = [0:pi/1000:pi/2];

ki = n_2*w/c;
kt = real(subs(N))*ki;

cos_phi = sqrt(1-sin(theta).^2/subs(N_squared));
r_TE = (cos(theta) - subs(N)*cos_phi)./(cos(theta) + subs(N)*cos_phi);
r_TM = (-subs(N) * cos(theta) + cos_phi)./(subs(N)*cos(theta) + cos_phi);

plot(theta, abs(r_TE).^2)
hold on
plot(theta, abs(r_TM).^2)

legend('Transverse Electric', 'Transverse Magnetic')
xlabel('Incidence Angle')
ylabel('Reflectance')

%yticks(1)

xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0' '\pi/8' '\pi/4' '3\pi/8' '\pi/2'})


%% Angle of Refraction
nexttile

% The frequencies of light for which to plot the angle of refraction

%ws_to_plot = [0.5*wp, wp, 1.5*wp];
ws_to_plot = [0.75*ws(1), ws(1), 1.25*ws(1)];

%legend('0.5\omega_p', '\omega_p', '1.5\omega_p')
legend('0.75\omega_1', '\omega_1', '1.25\omega_1')

for w = ws_to_plot
    subs(N_squared);
    
    ki = n_2*w/c;
    kt = real(subs(N))*ki;
    
    phi = asin(ki/kt * sin(theta));
    
    plot(theta, phi)
    hold on
end

xlabel('Incidence Angle')
ylabel('Transmission Angle')

xticks([0 pi/8 pi/4 3*pi/8 pi/2])
xticklabels({'0' '\pi/8' '\pi/4' '3\pi/8' '\pi/2'})

yticks([0 pi/8 pi/4 3*pi/8 pi/2])
yticklabels({'0' '\pi/8' '\pi/4' '3\pi/8' '\pi/2'})

%print -deps epsFig