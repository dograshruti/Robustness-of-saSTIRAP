clear t
clear cf
clear popk0
clear popk1
clear popk2
clear rmsarea
%%
tic
ntp =1*800;
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
%% Flags to be set by the user
deco=1;     % 0,1. 1- simulate with decoherence
sa=0;       % 0,1. 1- Turn on counterdiabatic term in the Hamiltonian
st=1;       % 0,1. 1- Turn on the STIRAP Hamiltonian
corr=0;     % 0,1. Always keep it 0
%% Parameter values to be defined by the user
sigma0 =13.571*10^(-9);
omega01 =44*10^6*2*pi;
omega12 =37*10^6*2*pi;
alpha=omega12/omega01;
w01=7.505*1e9*2*pi;
w12=7.234*1e9*2*pi;
%%%
mmf = 45;
sg = linspace(15, 45, mmf)*1e-9;
kkf = 61;
tsv = -linspace(1.5, 3, kkf);
%%%
wt=(w01+w12)/2;
phi01=0;
phi12=-0;
phi02=-pi/2-(phi01+phi12);
phit=-(pi+phi01+phi12+phi02)/2;
bigdelta=(w01-w12)/2;
%% Decay rates
GAMMA10=5e6;
GAMMA10phi=0e6;
GAMMA20phi=0e6;
GAMMA21=7e6;
GAMMA21phi=0e6;
gamma10=1/2*GAMMA10 + GAMMA10phi;
gamma20=1/2*GAMMA21 + GAMMA20phi;
gamma21=1/2*(GAMMA21 + GAMMA10) + GAMMA21phi;
%%
for mm=1:mmf
      sigma = sg(mm);
for kk=1:kkf
    
   rr = tsv(kk);
   ts=rr*sigma;
   bigdelta(mm,kk)=(w01-w12)/2;
   delta_phi01=linspace(-0.15,0.15,ntp)*(sg(mm)-rr*sg(mm))/4*1e8*0;
   delta_phi12=linspace(-0.2,0.2,ntp)*(sg(mm)-rr*sg(mm))/4*1e8*0;
   delta_phi01mm=[zeros(1,25),linspace(1,90,20)]*pi/180;
   delta_phi12mm=[zeros(1,25),linspace(1,90,20)]*pi/180;
%ti = -ni*sigma+ts;
%tf =nf*sigma;
ti = -3*sigma+ts;
tf = sigma*0.124;
%tf=ts/2;
nti(mm,kk)=ni;
ntf(mm,kk)=nf;
totaltime(mm,kk)=tf-ti;
dt = (tf - ti)/ntp;
delta01 = 0;
hb = 6.62607004*10^(-34)/2/pi;
% state
psi1 = [1; 0;0];
psi = psi1/norm(psi1);
rho0 = psi*psi';
    omega01_1=ones(1,25)*44*10^6*2*pi;
    omega01_2=(44+linspace(0,-2,20))*10^6*2*pi;
    omega01=[omega01_1, omega01_2];
    
    omega12_1=ones(1,25)*37*10^6*2*pi;
    omega12_2=(37+linspace(0,-2,20))*10^6*2*pi;
    omega12=[omega12_1, omega12_2];

%%%
for ii = 1:ntp
    rho_initial=rho0;
   
 t = ti + (ii - 1)*dt;
omega01t = omega01(mm)*exp(-(t)^2/2/sigma^2);
omega12t = omega12(mm)*exp(-(t - ts)^2/2/sigma^2);
omega02t = -2*ts*omega01t*omega12t/(omega01t^2+omega12t^2)/sigma^2;
tanq = omega01t/omega12t;
theta(ii,1) = atan(tanq);
%%%%
phi01t=phi01+2*sqrt(2)*atan(tanq)*corr+delta_phi01(ii)+delta_phi01mm(mm);
phi12t=phi12-5/sqrt(2)*atan(tanq)*corr+delta_phi12(ii)+delta_phi12mm(mm);
phi02t=phi02-(delta_phi01(ii)+delta_phi12(ii)+delta_phi01mm(mm)+delta_phi12mm(mm))+1/sqrt(2)*atan(tanq)*corr;
phit(mm,kk)=-(phi02t+pi)/2;
thetatotal(ii,2) = phi02t+phi01t+phi12t;
phi(mm,ii,1)=phi01t;
phi(mm,ii,2)=phi12t;
%%%
omega01tildet=sqrt(sqrt(2)*bigdelta(mm,kk)*omega02t);
o01(ii,1)=omega01t;
o12(ii,1)=omega12t;
o02(ii,1)=omega02t;
ot(ii,1)=omega01tildet;
%%%
H01=hb/2*[0, (omega01t)*exp(i*(phi01t)),0;
    (omega01t)*exp(-i*(phi01t)),delta01,0; 0,0,0];
%
H12=hb/2*[0,0,0; 0, 0, (omega12t)*exp(i*(phi12t));
    0, (omega12t)*exp(-i*(phi12t)),0];
%
H02=hb/2*[0,omega01tildet*exp(-i*bigdelta(mm,kk)*t+i*phit(mm,kk)),0;
    omega01tildet'*exp(i*bigdelta(mm,kk)*t-i*phit(mm,kk)),0,sqrt(2)*omega01tildet*exp(i*bigdelta(mm,kk)*t+i*phit(mm,kk));
    0,sqrt(2)*omega01tildet'*exp(-i*bigdelta(mm,kk)*t-i*phit(mm,kk)),0];
H=st*(H01+H12)+sa*H02;
dt0=dt;
k1=i/hb*(H*rho0 - rho0*H)*dt0;
%rho=rho0;
 %%
t = ti + (ii - 1)*dt + dt/2;
omega01t = omega01(mm)*exp(-(t)^2/2/sigma^2);
omega12t = omega12(mm)*exp(-(t - ts)^2/2/sigma^2);
omega02t = -2*ts*omega01t*omega12t/(omega01t^2+omega12t^2)/sigma^2;
tanq = omega01t/omega12t;
theta(ii,2) = atan(tanq);
%%%
phi01t=phi01+2*sqrt(2)*atan(tanq)*corr+delta_phi01(ii)+delta_phi01mm(mm);
phi12t=phi12-5/sqrt(2)*atan(tanq)*corr+delta_phi12(ii)+delta_phi12mm(mm);
phi02t=phi02-(delta_phi01(ii)+delta_phi12(ii)+delta_phi01mm(mm)+delta_phi12mm(mm))+1/sqrt(2)*atan(tanq)*corr;
phit(mm,kk)=-(phi02t+pi)/2;
thetatotal(ii,2) = phi02t+phi01t+phi12t;
Dphi(ii,1) = phi01t;
Dphi(ii,2) = phi12t;
Dphi(ii,3) = phi02t;
Dphi(ii,4) = phit(mm,kk);
omega01tildet=sqrt(sqrt(2)*bigdelta(mm,kk)*omega02t);
o01(ii,2)=omega01t;
o12(ii,2)=omega12t;
o02(ii,2)=omega02t;
ot(ii,2)=omega01tildet;
%%%
H01=hb/2*[0, (omega01t)*exp(i*(phi01t)),0;
    (omega01t)*exp(-i*(phi01t)),delta01,0; 0,0,0];
%
H12=hb/2*[0,0,0; 0, 0, (omega12t)*exp(i*(phi12t));
    0, (omega12t)*exp(-i*(phi12t)),0];
%
H02=hb/2*[0,omega01tildet*exp(-i*bigdelta(mm,kk)*t+i*phit(mm,kk)),0;
    omega01tildet'*exp(i*bigdelta(mm,kk)*t-i*phit(mm,kk)),0,sqrt(2)*omega01tildet*exp(i*bigdelta(mm,kk)*t+i*phit(mm,kk));
    0,sqrt(2)*omega01tildet'*exp(-i*bigdelta(mm,kk)*t-i*phit(mm,kk)),0];
H=st*(H01+H12)+sa*H02;
rho0 =  rho_initial + k1/2;
k2=i/hb*(H*rho0 - rho0*H)*dt0;
 %%
rho0 =  rho_initial + k2/2;
k3=i/hb*(H*rho0 - rho0*H)*dt0;
%rho=rho0;
 %%
t = ti + (ii - 1)*dt + dt;
omega01t = omega01(mm)*exp(-(t)^2/2/sigma^2);
omega12t = omega12(mm)*exp(-(t - ts)^2/2/sigma^2);
omega02t = -2*ts*omega01t*omega12t/(omega01t^2+omega12t^2)/sigma^2;
tanq = omega01t/omega12t;
theta(ii,3) = atan(tanq);
%%%
phi01t=phi01+2*sqrt(2)*atan(tanq)*corr+delta_phi01(ii)+delta_phi01mm(mm);
phi12t=phi12-5/sqrt(2)*atan(tanq)*corr+delta_phi12(ii)+delta_phi12mm(mm);
phi02t=phi02-(delta_phi01(ii)+delta_phi12(ii)+delta_phi01mm(mm)+delta_phi12mm(mm))+1/sqrt(2)*atan(tanq)*corr;
phit(mm,kk)=-(phi02t+pi)/2;
thetatotal(ii,1) = phi02t+phi01t+phi12t;
omega01tildet=sqrt(sqrt(2)*bigdelta(mm,kk)*omega02t);
o01(ii,3)=omega01t;
o12(ii,3)=omega12t;
o02(ii,3)=omega02t;
ot(ii,3)=omega01tildet;
%
H01=hb/2*[0, (omega01t)*exp(i*(phi01t)),0;
    (omega01t)*exp(-i*(phi01t)),delta01,0; 0,0,0];
%
H12=hb/2*[0,0,0; 0, 0, (omega12t)*exp(i*(phi12t));
    0, (omega12t)*exp(-i*(phi12t)),0];
%
H02=hb/2*[0,omega01tildet*exp(-i*bigdelta(mm,kk)*t+i*phit(mm,kk)),0;
    omega01tildet'*exp(i*bigdelta(mm,kk)*t-i*phit(mm,kk)),0,sqrt(2)*omega01tildet*exp(i*bigdelta(mm,kk)*t+i*phit(mm,kk));
    0,sqrt(2)*omega01tildet'*exp(-i*bigdelta(mm,kk)*t-i*phit(mm,kk)),0];
H=st*(H01+H12)+sa*H02;
rho0 =  rho_initial + k3;
%dt0=2*dt;
k4=i/hb*(H*rho0 - rho0*H)*dt0;
 %%
 rho1 =  rho_initial + 1/6*(k1 + 2*k2 + 2*k3 + k4);
 %%%
   Lrho=[GAMMA10*rho1(2,2), -gamma10*rho1(1,2), -gamma20*rho1(1,3);
     -gamma10*rho1(2,1), GAMMA21*rho1(3,3)-GAMMA10*rho1(2,2),-gamma21*rho1(2,3);
     -gamma20*rho1(3,1), -gamma21*rho1(3,2), -GAMMA21*rho1(3,3)];
  rho = rho1 + deco*Lrho*dt;
 %%%%%%%%%%%
 rho0 = rho;
 rhoii(:,:,ii)=rho;
 psiideal = [cos(theta(ii,1)); 0; -sin(theta(ii,1))];
 rhoideal(:,:,ii) = psiideal*psiideal';
 %rho(:,:,ii,mm,kk)=rho;
end
%%
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
end

popk0(mm,kk)=pop(ntp,1);
popk1(mm,kk)=pop(ntp,2);
popk2(mm,kk)=pop(ntp,3);
popk2max(mm,kk)=max(pop(:,3));
popk0optimal(mm,kk)=min(pop(:,1));
popk2optimal(mm,kk)=max(pop(:,3));
popk1optimal(mm,kk)=1-min(pop(:,1))-max(pop(:,3));
sensitivity(mm,kk)=max(pop(:,3))*(1-max(pop(:,3)))*2*rr/sigma*(2*0-rr);
%%%%%%%%%%%%%%%
end
%popm(mm,:)=pop(nt,:)
end
%%
fig3=figure 
[X,Y] = meshgrid(tsv,sg);
Z=popk2max;
surf(tsv, sg, Z);  shading flat;  %view([90 30]);
ylabel('\sigma (ns)');
xlabel('k');
title('State |2\rangle population');
ylim([sg(1) sg(end)]);
xlim([tsv(end) tsv(1)]);
zlim([0,1])
caxis([0 1])
colormap('jet')
colorbar;
%%
fig4=figure 
[X,Y] = meshgrid(tsv,sg);
Z=popk2;
%Z=totaltime;
h1=surf(X, Y, Z)
saveas(h1,sprintf('filename.jpg'))
%%
toc

%
