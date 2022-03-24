clc
clear t
clear cf
clear popk
tic
%%
nt =1*1801;
ntp=ceil(nt)-1;

% Flags
st=1;         % 1- turns on STIRAP Hamiltonian
sa=1;         % 1- turn on the counterdiabatic part of the Hamiltonian
deco=1;        % 1- sets in decoherence
direct02coupling=0; % 1 for directly coupled 0-2 drive
corr=0;        % keep at 0
area_plot=1;   % To be kept 1 for area map
%%%%%%%%%%%%%%%
theta=zeros(ntp,3);
theta_av=zeros(1,ntp);
o01=zeros(ntp,3);
o12=zeros(ntp,3);
o02=zeros(ntp,3);
ot=zeros(ntp,3);
rhoii=zeros(3,3,ntp);
rhoideal=zeros(3,3,ntp);
pop=zeros(ntp,3);
popideal=zeros(ntp,3);
dq=zeros(ntp,1);
sq=zeros(ntp,2);
purity=zeros(ntp,1);
%% Example 1
omega01 = 44*10^6*2*pi;%44*10^6*2*pi;%46.3*10^6*2*pi;%41*10^6*2*pi;
omega12 = 37*10^6*2*pi;%37*10^6*2*pi;%46.3*10^6*2*pi;%35*10^6*2*pi;
omega02 = 40.7416*10^6*2*pi; %47*10^6*2*pi;
alpha=omega12/omega01;

%%%
w01=7.394*10^(9)*2*pi;
w12=7.099*10^(9)*2*pi;
wt=(w01+w12)/2;
phi01=0;
phi12=-0;
phi02=-pi/2-(phi01+phi12);
phit=-(pi+phi01+phi12+phi02)/2;
%phit=-pi/4;
bigdelta=(w01-w12)/2;
% from robustness paper
bigdelta = 293*1e6/2*2*pi; 
%% %%%%%%%%% User defined parameters %%%%%%%%%%%%
mmf=31;
% sigma_range=linspace(sigma0,40*1e-9,mmf);
kkf=41; 
% ts_range=-linspace(3,1.5,kkf);
sigma=30e-9;
ts=-1.5*sigma;  
omega02_pi=-ts/sigma/sigma;
omega01_range=2*linspace(0.0,40,kkf)*2*pi*1e6;
omega12_range=2*linspace(0.0,29,kkf)*2*pi*1e6; 
omega02_range=1*linspace(0,sqrt(6)*omega02_pi,mmf);

for mm=1:mmf
    
for kk=1:kkf
      if(area_plot==1)
        omega01 = omega01_range(kk); %0 - 40*10^6*2*pi;
        omega12 = omega12_range(kk); %0 - 26*10^6*2*pi;
        omega02 = omega02_range(mm); %0 - 94*10^6*2*pi;
        alpha=omega12/omega01;
    end
        rr=ts/sigma;
    ni=3;
    nf=3;
ti = -ni*sigma+ts;
tf =nf*sigma;
dt = (tf - ti)/(nt-1);
tt(kk,1)=tf-ti;
delta01 = 0;
hb = 6.62607004*10^(-34)/2/pi;

ti01=-ni*sigma;
tf01= nf*sigma;
ti12=-ni*sigma+ts;
tf12= nf*sigma+ts;
%%
    %%%% Evaluate area %%%%%%
fun_area = @(q) sqrt(omega01^2*exp(-(q).^2/sigma^2) + omega12^2*exp(-(q-ts).^2/sigma^2));

areaA = integral(fun_area,ti,tf);
Arms(kk,1)=areaA;

area_range=linspace(0,6,mmf);
%

fun_area02 = @(q) (omega02/omega02_pi)^2./cosh(-ts/sigma^2*(q-ts/2))*(-ts/sigma^2)/sqrt(2);

if(direct02coupling==1)
fun_area02 = @(q) -area_range(mm)*2*ts*omega01*omega12*exp(-(q).^2/2/sigma^2).*exp(-(q - ts).^2/2/sigma^2)./(omega01^2*exp(-(q).^2/sigma^2)+omega12^2*exp(-(q - ts).^2/sigma^2))/sigma^2;
end

areaA02=integral(fun_area02,ti,tf);
A02(mm,1) = areaA02;

% A02(mm,1) = areaA02/(sqrt(2)*bigdelta)*(-ts)/sigma^2;
%%
% Decay rates
GAMMA10=5e6;
GAMMA10phi=0e6;
GAMMA20phi=0e6;
GAMMA21=7e6;
GAMMA21phi=0e6;
gamma10=1/2*GAMMA10 + GAMMA10phi;
gamma20=1/2*GAMMA21 + GAMMA20phi;
gamma21=1/2*(GAMMA21 + GAMMA10) + GAMMA21phi;
% state
psi1 = [1; 0;0];
psi = psi1/norm(psi1);

rho0 = psi*psi';
for ii = 1:ntp 
    rho_initial=rho0;
   
 t = ti + (ii - 1)*dt;
 
omega01t = omega01*exp(-(t)^2/2/sigma^2);
omega12t = omega12*exp(-(t - ts)^2/2/sigma^2);
omega02t = -2*ts*omega01t*omega12t/(omega01t^2+omega12t^2)/sigma^2*area_range(mm);
omega02exp=-ts/sigma^2/cosh(-ts/sigma^2*(t-ts/2))*area_range(mm);
omega2pht=sqrt(sqrt(2)*bigdelta*omega02exp);
tanq = omega01t/omega12t;
theta(ii,1) = atan(tanq);
%%%%%%%%%%%%
if(t>=ti01 && t<=tf01+dt)
omega01t = omega01t;
 else
     omega01t=0;
 end
 if(t>=ti12 && t<=tf12+dt)
omega12t = omega12t;
else
     omega12t=0;
 end
%%%%
omega01tildet=sqrt(sqrt(2)*bigdelta*omega02t);
o01(ii,1)=omega01t;
o12(ii,1)=omega12t;
o02(ii,1)=omega02t;
ot(ii,1)=omega01tildet;
o2pht(ii,1)=omega2pht;
%%%
H01=hb/2*[0, (omega01t)*exp(i*phi01),0;
    (omega01t)*exp(-i*phi01),delta01,0; 0,0,0];
H12=hb/2*[0,0,0; 0, 0, (omega12t)*exp(i*phi12);
    0, (omega12t)*exp(-i*phi12),0];
H02=hb/2*[0,omega2pht*exp(-i*bigdelta*t+i*phit),0;
     omega2pht'*exp(i*bigdelta*t-i*phit),0,sqrt(2)*omega2pht*exp(i*bigdelta*t+i*phit);
     0,sqrt(2)*omega2pht'*exp(-i*bigdelta*t-i*phit),0];
 
 if(direct02coupling==1)
     H02=hb/2*[0,0,omega02t*exp(i*phi02); 0,0,0; omega02t*exp(-i*phi02),0,0];
 end

H=st*(H01+H12)+sa*H02;
dt0=dt;
k1=i/hb*(H*rho0 - rho0*H)*dt0;
 %%
t = ti + (ii - 1)*dt + dt/2;

omega01t = omega01*exp(-(t)^2/2/sigma^2);
omega12t = omega12*exp(-(t - ts)^2/2/sigma^2);
omega02t = -2*ts*omega01t*omega12t/(omega01t^2+omega12t^2)/sigma^2*area_range(mm);
omega02exp=-ts/sigma^2*1/cosh(-ts/sigma^2*(t-ts/2))*area_range(mm);
omega2pht=sqrt(sqrt(2)*bigdelta*omega02exp);
tanq = omega01t/omega12t;
theta(ii,2) = atan(tanq);
%%%
 if(t>=ti01 && t<=tf01+dt)
omega01t = omega01t;
 else
     omega01t=0;
 end
 if(t>=ti12 && t<=tf12+dt)
omega12t = omega12t;
else
     omega12t=0;
 end
omega01tildet=sqrt(sqrt(2)*bigdelta*omega02t);
o01(ii,2)=omega01t;
o12(ii,2)=omega12t;
o02(ii,2)=omega02t;
ot(ii,2)=omega01tildet;
o2pht(ii,2)=omega2pht;
%%%
H01=hb/2*[0, (omega01t)*exp(i*phi01),0;
    (omega01t)*exp(-i*phi01),delta01,0; 0,0,0];
H12=hb/2*[0,0,0; 0, 0, (omega12t)*exp(i*phi12);
    0, (omega12t)*exp(-i*phi12),0];
H02=hb/2*[0,omega2pht*exp(-i*bigdelta*t+i*phit),0;
     omega2pht'*exp(i*bigdelta*t-i*phit),0,sqrt(2)*omega2pht*exp(i*bigdelta*t+i*phit);
     0,sqrt(2)*omega2pht'*exp(-i*bigdelta*t-i*phit),0];
 
  if(direct02coupling==1)
     H02=hb/2*[0,0,omega02t*exp(i*phi02); 0,0,0; omega02t*exp(-i*phi02),0,0];
 end

H=st*(H01+H12)+sa*H02;
rho0 =  rho_initial + k1/2;
k2=i/hb*(H*rho0 - rho0*H)*dt0;
 %%
rho0 =  rho_initial + k2/2;
k3=i/hb*(H*rho0 - rho0*H)*dt0;
 %%
t = ti + (ii - 1)*dt + dt;

omega01t = omega01*exp(-(t)^2/2/sigma^2);
omega12t = omega12*exp(-(t - ts)^2/2/sigma^2);
omega02t = -2*ts*omega01t*omega12t/(omega01t^2+omega12t^2)/sigma^2*area_range(mm);
omega02exp=-ts/sigma^2*1/cosh(-ts/sigma^2*(t-ts/2))*area_range(mm);
omega2pht=sqrt(sqrt(2)*bigdelta*omega02exp);
tanq = omega01t/omega12t;
theta(ii,3) = atan(tanq);
%%%
 if(t>=ti01 && t<=tf01+dt)
omega01t = omega01t;
 else
     omega01t=0;
 end
 if(t>=ti12 && t<=tf12+dt)
omega12t = omega12t;
else
     omega12t=0;
 end
omega01tildet=sqrt(sqrt(2)*bigdelta*omega02t);
o01(ii,3)=omega01t;
o12(ii,3)=omega12t;
o02(ii,3)=omega02t;
ot(ii,3)=omega01tildet;
%
H01=hb/2*[0, (omega01t)*exp(i*phi01),0;
    (omega01t)*exp(-i*phi01),delta01,0; 0,0,0];
H12=hb/2*[0,0,0; 0, 0, (omega12t)*exp(i*phi12);
    0, (omega12t)*exp(-1i*phi12),0];
H02=hb/2*[0,omega2pht*exp(-i*bigdelta*t+i*phit),0;
     omega2pht'*exp(i*bigdelta*t-i*phit),0,sqrt(2)*omega2pht*exp(i*bigdelta*t+i*phit);
     0,sqrt(2)*omega2pht'*exp(-i*bigdelta*t-i*phit),0];
 
  if(direct02coupling==1)
     H02=hb/2*[0,0,omega02t*exp(i*phi02); 0,0,0; omega02t*exp(-i*phi02),0,0];
 end

H=st*(H01+H12)+sa*H02;
rho0 =  rho_initial + k3;
k4=i/hb*(H*rho0 - rho0*H)*dt0;
 %%
 rho1 =  rho_initial + 1/6*(k1 + 2*k2 + 2*k3 + k4);
 
 Lrho=[GAMMA10*rho1(2,2), -gamma10*rho1(1,2), -gamma20*rho1(1,3);
     -gamma10*rho1(2,1), GAMMA21*rho1(3,3)-GAMMA10*rho1(2,2),-gamma21*rho1(2,3);
     -gamma20*rho1(3,1), -gamma21*rho1(3,2), -GAMMA21*rho1(3,3)];
  rho = rho1 + deco*Lrho*dt;
 
 rho0 = rho;
 rhoii(:,:,ii)=rho;
 rho_iimmkk(:,:,ii,mm,kk)=rho;
 psiideal = [cos(theta(ii,1)); 0; -sin(theta(ii,1))];
 rhoideal(:,:,ii) = psiideal*psiideal';
end
for jj=1:ntp
    pop(jj,1)=real(rhoii(1,1,jj));
    pop(jj,2)=real(rhoii(2,2,jj));
    pop(jj,3)=real(rhoii(3,3,jj));
    dq(jj,1)=rhoii(1,3,jj);
    sq(jj,1)=rhoii(1,2,jj);
    sq(jj,2)=rhoii(2,3,jj);
    popideal(jj,1)=real(rhoideal(1,1,jj));
    popideal(jj,2)=real(rhoideal(2,2,jj));
    popideal(jj,3)=real(rhoideal(3,3,jj));
    purity(jj,1)=trace(rhoii(:,:,jj)*rhoii(:,:,jj));
%
end
pp=0.57;
popk(mm,kk)=pop(ntp,3);
popkmax(mm,kk)=max(pop(:,3));
popkchop1(mm,kk)=pop(ceil(pp*ntp),3);
end
%
end
%%
pp=1200/1800
if (area_plot==1)
    
figure

X1=area_range; 
Y1=Arms/pi; %+[0.1; 0.2];
Z=abs(squeeze(rho_iimmkk(3,3,ceil(pp*ntp),:,:))); %popk;
[X,Y] = meshgrid(Y1,X1);
h1=surf(X, Y, Z)
view([90 -90]);
caxis([0,1])
colormap('hot')
colorbar
shading 'flat'
xlabel('A');
ylabel('A_{02}/\pi')
end
%%
 toc

