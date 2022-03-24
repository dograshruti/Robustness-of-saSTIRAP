function F = evolut(t, x, st, sa, dc02,params)
syms t1 real;

% params=[sigma; ts; omega01; omega12; ti; tf; angl; det; phi];

sigma = params(1); 
ts = params(2); 
omega01 = params(3); 
omega12 = params(4); 
ti = params(5); 
tf = params(6); 
angl = params(7); 
det = params(8); 
dynamic_phase = params(9); 
cosh02=params(10);
phi01=params(11);
phi12=params(12);
phi20=params(13);
phi_kk=params(14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega01t = omega01*exp(-(t1)^2/2/sigma^2);
omega12t = omega12*exp(-(t1-ts)^2/2/sigma^2);
omega02t = -2*ts/sigma^2*omega01t*omega12t/(omega01t^2+omega12t^2)*angl;
if(cosh02==1)
omega02t=-ts/sigma^2/cosh(-ts/sigma^2*(t1-ts/2))*angl;
end
omega2pht = sqrt(sqrt(2)*det*omega02t);%*angl;

hb = 1;% 6.62607004*10^(-34)/2/pi;
Thetat = atan(omega01t/omega12t)*sqrt(angl);
phi01t = phi01 + (dynamic_phase)*4/sqrt(2)*hb*Thetat;
phi12t = phi12 - (dynamic_phase)*5/sqrt(2)*hb*Thetat;
phi20t = phi20 + (dynamic_phase)*hb/sqrt(2)*Thetat;
phi2pht = 1*(-(phi20t+pi)/2 + phi_kk);
%phi2pht = -(phi20t+pi)/2; 

H01=hb/2*[0, omega01t*exp(i*phi01t),0; omega01t*exp(-i*phi01t),0,0; 0,0,0];
H12=hb/2*[0,0,0; 0,0,omega12t*exp(i*phi12t); 0,omega12t*exp(-i*phi12t),0];
H02=hb/2*[0,1*omega2pht*exp(i*(phi2pht-det*t1)),0;...
   1*omega2pht*exp(-i*(phi2pht-det*t1)),0*det,sqrt(2)*omega2pht*exp(i*(phi2pht+det*t1));...
   0,sqrt(2)*omega2pht*exp(-i*(phi2pht+det*t1)),0*det];

if dc02==1
    H02=hb/2*[0,0,exp(-i*phi20t)*omega02t; 0,0,0; exp(i*phi20t)*omega02t,0,0];
end
    
H=st*(H01+H12)+sa*H02;
H_sai = subs(H,t1,t);
F = double(-1i*H_sai*x);
%%%%%%%%%%%%%%

end
