%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user parameters have same names as in the article 
%"Experimental demonstration of robustness under scaling errors for 
% superadiabatic population transfer in a superconducting circuit, 
% Shruti Dogra, Antti Vepsäläinen, and Gheorghe Sorin Paraoanu".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It is important to set the correct flags for a desired simulation
% Description of each flag is given separately.
% Flags can assume values 0 or 1. General rule:
% 1- Feature ON (descriptions are given below), 0- Feature OFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega01=2*pi*8e6; 
omega12=2*pi*8e6; 
sigma=30e-9;
ts=-1.5*sigma;
ti=-3*sigma+ts;
tf=3*sigma;
angl=1;     % multiple of pi, such that Area A_{02}=angl*pi
v01=7.2e9;
v12=6.9e9;
det=2*pi*(v01-v12)/2;
phi_range=0;
phi01=0;
phi12=0;
phi20=-pi/2;
interval = [ti tf];

%%%% Flags
st=1        % 1,0. 1- set STIRAP Hamiltonian
sa=1        % 1,0. 1- set Counterdiabatic drive Hamiltonian
dc02=0      % 1,0. 1- direct coupling between levels 0-2
cosh02=0    % 1,0. 1- omega_{02}(t) is approximated to cosh pulse
area_plot=1 % 1,0. 1- set STIRAP/saSTIRAP p_2 surface map vs A and A_{02}
rabi02=0    % 1,0. 1- Acquire Rabi oscillations in the 0-2 subspace
dynamic_phase=0 % 1,0. Always keep it 0
operator_contruction=0; % 1,0. 1- Constructs the unitary operator generated by the chosen Hamiltonian
%%%%%%%%%%
tic


for kk=1:length(phi_range)
phi_kk = phi_range(kk);

%%%%%%
if (rabi02==1)
    angl_range=linspace(0,10,101);
else angl_range=angl;
end
%%%%%%%%%%
if(area_plot==1)
    angl_range=linspace(0,8,41);
    omega_range=2*pi*linspace(0.001,30,31)*1e6;
    omega12_range=2*pi*linspace(0.001,30,31)*1e6;
elseif (area_plot==0 && rabi02==0)
    angl_range=angl;
    omega_range=omega01;
elseif (area_plot==0 && rabi02==1)
    omega_range=omega01;
end
%%%%%%%%%%%
for jj=1:length(omega_range)
    
    if (area_plot==1)
    omega01=omega_range(jj);
    omega12=omega12_range(jj);
    end
    
for ii=1:length(angl_range)
    
        angl=angl_range(ii);

params=[sigma; ts; omega01; omega12; ti; tf; angl; det; dynamic_phase; cosh02; phi01; phi12; phi20; phi_kk];
%
%options_dif = odeset('AbsTol', 3.8e-9,'RelTol', 3.8e-7,'Refine', 10);
options_dif = odeset('AbsTol', 0.5e-5,'RelTol', 1e-5,'Refine', 1);

% Schrodinger's equation for pure state evolution

x1 = [1; 0; 0]; % initial state vector
[t x] = ode45(@(t,x) evolut(t,x,st,sa,dc02,params), interval, x1,options_dif);

xr=x.*conj(x);
population1(ii,jj,kk)=xr(end,1);
population2(ii,jj,kk)=xr(end,2);
population3(ii,jj,kk)=xr(end,3);

state1(ii,jj,kk)=x(end,1);
state2(ii,jj,kk)=x(end,2);
state3(ii,jj,kk)=x(end,3);

if(operator_contruction==1)
    [t2 x2] = ode45(@(t,x) evolut(t,x,st,sa,dc02,params), interval, [0;1;0],options_dif);
    [t3 x3] = ode45(@(t,x) evolut(t,x,st,sa,dc02,params), interval, [0;0;1],options_dif);
    %for ijk=1:length(t)
    xn = [x(end,:);x2(end,:);x3(end,:)]
    % end
end
%ii
end
jj;
end
kk;
end
toc
%%
 figure()
 hold on
if(rabi02==1)   
plot(angl_range, population1);
hold on
plot(angl_range, population2);
hold on
plot(angl_range, population3);
elseif (area_plot==1)
surf(squeeze(population3(:,:)))
shading 'flat'
else
    plot(t*1e9,xr)
end
%%
if rabi02==0
figure
tiledlayout(3,1)
nexttile
plot(t*1e9,xr)
    xlabel('Time [ns]');
    ylabel('Populations');
    legend('p_0','p_1','p_2');
    title("st=" + st + ", sa=" + sa + ", dc02=" + dc02 + ", A_{02}=" + angl + "\pi" ...
        + ", \Omega_{01}=" + omega01*1e-6/2/pi + "MHz" ...
        + ", \Omega_{12}=" + omega12*1e-6/2/pi + "MHz"...
        + ", \sigma=" + sigma*1e9 + "ns" + ", t_s=" + ts*1e9 + "ns")

nexttile
plot(t*1e9,real(x))
    xlabel('Time [ns]');
    ylabel('Coefficients (Real)');
    legend('c_0','c_1','c_2');
    
nexttile
plot(t*1e9,imag(x))
    xlabel('Time [ns]');
    ylabel('Coefficients (Imag)');
    legend('c_0','c_1','c_2');
    
else
figure
tiledlayout(1,2)
nexttile
plot(t*1e9,real(x))
nexttile
plot(t*1e9,imag(x))
end
% 
%  %%
% figure
% tiledlayout(3,1)
% 
% nexttile
% X1=angl_range.^1;
% Y1=omega_range/2/pi*1e-6;
% Z=squeeze(population3(:,:));
% [X,Y] = meshgrid(Y1,X1);
% h1=surf(X, Y, Z)
% view([90 -90]);
% caxis([0,1])
% colormap('hot')
% colorbar
% shading 'flat'
% %
% 
% nexttile
% %X1=angl_range;
% Y1=omega_range/2/pi*1e-6;
% Z=squeeze(population2(:,:));
% [X,Y] = meshgrid(Y1,X1);
% h1=surf(X, Y, Z)
% view([90 -90]);
% caxis([0,1])
% colormap('hot')
% colorbar
% shading 'flat'
% 
% 
% nexttile
% %X1=angl_range;
% Y1=omega_range/2/pi*1e-6;
% Z=squeeze(population3(:,:));
% [X,Y] = meshgrid(Y1,X1);
% h1=surf(X, Y, Z)
% view([90 -90]);
% caxis([0,1])
% colormap('hot')
% colorbar
% shading 'flat'

