%%-------------------------------------------------------------------
% this code is for calculating physical properties and matrices for 
% the simulation & experimental setup
% by Yashwanth Nakka
% all dimensions are in SI
clc;
tic;
digits(5);
%%% for simulation 
%%Beam Properties
w=3.6*10^-2; %width of the beam and the piezo 
rho_b=2758; %density of the beam & vicon markers
lb=29.7*10^-2; % length of the beam
hb=0.45*10^-3; % thickness of the beam
Eb= 7.0e+10;% yound modulus of the beam

%% Piezo properties
rho_p=7800;% density of piezo
lp=7.24*10^-2; % length of the piezo
hp=0.48*10^-3 ; % thickness of the piezo
Ep=52*10^9;% yound modulus of the piezo
d31 =-190*10^-12; % peizo
% Piezo1 location 
lp_1=1.1*10^-2; % distance of the root of the piezo from that beamroot
lp_2=lp_1+lp;
% distance between Piezo 1 and Piezo 2 
dpp=0.47*10^-2;
lp1_1= lp_2 +dpp;
lp1_2= lp1_1 +lp;

%% Ball Bearing,Mass Moment of Inertia of Inner race Properties
In_in= 5*10^-3;% inner radius
In_out= 8*10^-3;% outer radius
rho_ball=7990 ;% density of steel 
w_ball=8*10^-3; % thickness of the ball bearing 
mass_ball= rho_ball*pi*(In_out^2 - In_in^2)*w_ball;% Mass of Ball bearing inner race
Inertia_ball=mass_ball*(In_out^2 + In_in^2)*0.5; % Mass Inertial of the ball bearing

%% Central Hub properties
cyl_in= 0;% inner radius
cyl_out= 5*10^-3;% outer radius
rho_cyl=2700 ;% density of alumininum 
cyl_length=10*10^-2; % length of the cylinder
mass_cyl= (rho_cyl*pi*(cyl_out^2 - cyl_in^2)*cyl_length) ; % Mass of cylinder 
Inertia_cyl=mass_cyl*(cyl_out^2 + cyl_in^2)*0.5; % Mass Inertial of the cylinder

%% Central Hub Flange properties
Flange_in= 5*10^-3;% inner radius
Flange_out= 6*10^-3;% outer radius
rho_Flange=2700 ;% density of alumininum 
Flange_length=4*10^-2; % length of the Flange projection
mass_Flange= rho_Flange*pi*(Flange_out^2 - Flange_in^2)*Flange_length; % Mass of Flange 
Inertia_Flange=mass_Flange*(Flange_out^2 + Flange_in^2)*0.5; % Mass Inertial of the Flange

%% vicon markers

mv= 2*0.40*10^-3;
lv= [5*0.01 10*0.01 15*0.01 20*0.01 25*0.01 30*0.01];


%% Total Mass Inertia of the hub includes cylinder and one bearing
% need to add one more bearing
Total_J = 2*Inertia_ball +Inertia_cyl + Inertia_Flange;% +((0.65*10^-3)*(3*10^-3)^2)/2;




%% properties at the location of piezo 
% neutral axis at the location of the piezo 
Hn= ((Ep*(hp^2)/2) +(Eb*hb*((hb/2) + hp)))/(Ep*hp + Eb*hb);
% total rigidity of the beam for the new neutral axis
EIb = Eb*((w*hb^3)/12) ;
EIp =  Eb*(w*hb *(hp + (0.5*hb) -Hn)^2)+ Ep*((w*hp^3)/12) + Ep*(w*hp *(Hn - (hp/2))^2);

% the total rigidity
EI= EIp + EIb;

% total mass per unit length
mr_beam= 0.040426488;%rho_b*hb*w;
mr_total= rho_p*hp*w + 0.040426488;% rho_b*hb*w;

%%modes and derivatives computation for 
syms 'x';
modes=4;% number of modes
u= sym(zeros(modes,1));
du= sym(zeros(modes,1));
ddu= sym(zeros(modes,1));
ddddu= sym(zeros(modes,1));
for k=1:modes
    cons= k*(1+k)*(2+k)*(3+k);
    %var=x/lb;
    %u(k)= (var^(1+k))*(6+(k^2*(1-var)^2)+ (k*(5 - 6*var + var^2)));
   u(k)=1-(cos(k*pi*x/lb))+(0.5*((-1)^(k+1))*((k*pi*x/lb)^2));
     %u(k)= (x/lb)^(k+1);
    du(k)=diff(u(k),x);
    ddu(k)=diff(diff(u(k),x),x);
    ddddu(k)= diff(diff(diff(diff(u(k),x),x),x),x); 
end

%% computation of system matrices
% 
% piecewise functions for mass per unit length and EI
syms 's';
p=s;
% p1= mr_beam*int(p,s,Flange_out+x,Flange_out+lb);
% p2= (mr_total*int(p,s,Flange_out+x,Flange_out+lp1_2))...
%     + (mr_beam*int(p,s,Flange_out+lp1_2,Flange_out+lb));
% p3= (mr_beam*int(p,s,Flange_out+x,Flange_out+lp1_1))...
%     +(mr_total*int(p,s,Flange_out+lp1_1,Flange_out+lp1_2))...
%     + (mr_beam*int(p,s,Flange_out+lp1_2,Flange_out+lb));
% p4= (mr_total*int(p,s,Flange_out+x,Flange_out+lp_2))...
%     +(mr_beam*int(p,s,Flange_out+lp_2,Flange_out+lp1_1))...
%     +(mr_total*int(p,s,Flange_out+lp1_1,Flange_out+lp1_2))...
%     + (mr_beam*int(p,s,Flange_out+lp1_2,Flange_out+lb));
% p5= (mr_beam*int(p,s,Flange_out+x,Flange_out+lp_1))...
%     +(mr_total*int(p,s,Flange_out+lp_1,Flange_out+lp_2))...
%     +(mr_beam*int(p,s,Flange_out+lp_2,Flange_out+lp1_1))...
%     +(mr_total*int(p,s,Flange_out+lp1_1,Flange_out+lp1_2))...
%     + (mr_beam*int(p,s,Flange_out+lp1_2,Flange_out+lb));


p10=  mr_beam*int(p,s,Flange_out+x,Flange_out+lb) +mv*(Flange_out+lb);

p9= mr_beam*int(p,s,Flange_out+x,Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);

p8=  +mr_beam*int(p,s,Flange_out +x,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);

p7=  mr_total*int(p,s,Flange_out+x,Flange_out+lp1_2)...
    +mr_beam*int(p,s,Flange_out +lp1_2,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);

p6=  mr_total*int(p,s,Flange_out+x,Flange_out+lv(3)) +mv*(Flange_out+lv(3)) ...  
    +mr_total*int(p,s,Flange_out+lv(3),Flange_out+lp1_2)...
    +mr_beam*int(p,s,Flange_out +lp1_2,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);

p5= mr_total*int(p,s,Flange_out+x,Flange_out+lv(2)) + mv*(Flange_out+lv(2)) ...
    +mr_total*int(p,s,Flange_out+lv(2),Flange_out+lv(3)) +mv*(Flange_out+lv(3)) ...  
    +mr_total*int(p,s,Flange_out+lv(3),Flange_out+lp1_2)...
    +mr_beam*int(p,s,Flange_out +lp1_2,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);

p4= mr_beam*int(p,s,Flange_out+x,Flange_out+lp1_1)...
    +mr_total*int(p,s,Flange_out+lp1_1,Flange_out+lv(2)) + mv*(Flange_out+lv(2)) ...
    +mr_total*int(p,s,Flange_out+lv(2),Flange_out+lv(3)) +mv*(Flange_out+lv(3)) ...  
    +mr_total*int(p,s,Flange_out+lv(3),Flange_out+lp1_2)...
    +mr_beam*int(p,s,Flange_out +lp1_2,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);

p3= mr_total*int(p,s,Flange_out+x,Flange_out+lp_2)...
    +mr_beam*int(p,s,Flange_out+lp_2,Flange_out+lp1_1)...
    +mr_total*int(p,s,Flange_out+lp1_1,Flange_out+lv(2)) + mv*(Flange_out+lv(2)) ...
    +mr_total*int(p,s,Flange_out+lv(2),Flange_out+lv(3)) +mv*(Flange_out+lv(3)) ...  
    +mr_total*int(p,s,Flange_out+lv(3),Flange_out+lp1_2)...
    +mr_beam*int(p,s,Flange_out +lp1_2,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);
 

p2= mr_total*int(p,s,Flange_out+x,Flange_out+lv(1)) + mv*(Flange_out+lv(1)) ...
    +mr_total*int(p,s,Flange_out+lv(1),Flange_out+lp_2)...
    +mr_beam*int(p,s,Flange_out+lp_2,Flange_out+lp1_1)...
    +mr_total*int(p,s,Flange_out+lp1_1,Flange_out+lv(2)) + mv*(Flange_out+lv(2)) ...
    +mr_total*int(p,s,Flange_out+lv(2),Flange_out+lv(3)) +mv*(Flange_out+lv(3)) ...  
    +mr_total*int(p,s,Flange_out+lv(3),Flange_out+lp1_2)...
    +mr_beam*int(p,s,Flange_out +lp1_2,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);

p1= mr_beam*int(p,s,Flange_out+x,Flange_out+lp_1)...
    +mr_total*int(p,s,Flange_out+lp_1,Flange_out+lv(1)) + mv*(Flange_out+lv(1)) ...
    +mr_total*int(p,s,Flange_out+lv(1),Flange_out+lp_2)...
    +mr_beam*int(p,s,Flange_out+lp_2,Flange_out+lp1_1)...
    +mr_total*int(p,s,Flange_out+lp1_1,Flange_out+lv(2)) + mv*(Flange_out+lv(2)) ...
    +mr_total*int(p,s,Flange_out+lv(2),Flange_out+lv(3)) +mv*(Flange_out+lv(3)) ...  
    +mr_total*int(p,s,Flange_out+lv(3),Flange_out+lp1_2)...
    +mr_beam*int(p,s,Flange_out +lp1_2,Flange_out+lv(4)) + mv*(Flange_out+lv(4))...
    +mr_beam*int(p,s,Flange_out+lv(4),Flange_out+lv(5)) +mv*(Flange_out+lv(5))...
    +mr_beam*int(p,s,Flange_out+lv(5),Flange_out+lb) +mv*(Flange_out+lb);



% integrations required for matrices
A=zeros(modes,modes);
B=zeros(modes,modes);
C=zeros(1,modes);
D=zeros(modes,modes);
E=zeros(modes,modes);
Q=0;
for i=1:modes
for j=1:modes
phi=u(i)*u(j);
phip=du(i)*du(j);
phie=u(i)*ddddu(j);
A1= int(phi,x,0,lp_1);
A2= int(phi,x,lp_1,lp_2);
A3= int(phi,x,lp_2,lp1_1);
A4= int(phi,x,lp1_1,lp1_2);
A5= int(phi,x,lp1_2,lb);

A(i,j)= (mr_beam*A1) + (mr_total*A2) + (mr_beam*A3) + (mr_total*A4) + (mr_beam*A5)...
    + mv*subs(phi,x,lv(1))+ mv*subs(phi,x,lv(2))+ mv*subs(phi,x,lv(3))+ mv*subs(phi,x,lv(4))...
    + mv*subs(phi,x,lv(5))+ mv*subs(phi,x,lv(6));

% B(i,j)= int(p5*phip,x,0,lp_1) + int(p4*phip,x,lp_1,lp_2)...
%     + int(p3*phip,x,lp_2,lp1_1) + int(p2*phip,x,lp1_1,lp1_2) + int(p1*phip,x,lp1_2,lb);

B(i,j) = int(p1*phip,x,0,lp_1) +int(p2*phip,x,lp_1,lv(1))...
    +int(p3*phip,x,lv(1),lp_2) +int(p4*phip,x,lp_2,lp1_1)...
    +int(p5*phip,x,lp1_1,lv(2)) + int(p6*phip,x,lv(2),lv(3))...
    +int(p7*phip,x,lv(3),lp1_2) +int(p8*phip,x,lp1_2,lv(4))...
    +int(p9*phip,x,lv(4),lv(5)) +int(p10*phip,x,lv(5),lv(6));


E(i,j)= int(EIb*phie,x,0,lp_1) + int(EI*phie,x,lp_1,lp_2)...
    + int(EIb*phie,x,lp_2,lp1_1)+int(EI*phie,x,lp1_1,lp1_2)+int(EIb*phie,x,lp1_2,lb);

D(i,j)= A1+A2+A3+A4+A5; 
end

C(i)= int(mr_beam*(x+Flange_out)*u(i),x,0,lp_1)...
    +int(mr_total*(x+Flange_out)*u(i),x,lp_1,lp_2)...
    +int(mr_beam*(x+Flange_out)*u(i),x,lp_2,lp1_1)...
    +int(mr_total*(x+Flange_out)*u(i),x,lp1_1,lp1_2)...
    +int(mr_beam*(x+Flange_out)*u(i),x,lp1_2,lb)...
   + mv*subs(u(i)*(x+Flange_out),x,lv(1))+ mv*subs(u(i)*(x+Flange_out),x,lv(2))...
    + mv*subs(u(i)*(x+Flange_out),x,lv(3))+ mv*subs(u(i)*(x+Flange_out),x,lv(4))...
    + mv*subs(u(i)*(x+Flange_out),x,lv(5))+ mv*subs(u(i)*(x+Flange_out),x,lv(6));
end

% inertia due to vicon markers
vicon= (0.40*10^-3)*((5+0.6)^2+(10+0.6)^2+(15+0.6)^2+(20+0.6)^2+(25+0.6)^2+(30+0.6)^2)*10^-4;
Q= Total_J ...
    + 2*vpa(int(mr_beam*(x+Flange_out)^2,x,0,lp_1)...
    +int(mr_total*(x+Flange_out)^2,x,lp_1,lp_2)...
    +int(mr_beam*(x+Flange_out)^2,x,lp_2,lp1_1)...
    +int(mr_total*(x+Flange_out)^2,x,lp1_1,lp1_2)...
    +int(mr_beam*(x+Flange_out)^2,x,lp1_2,lb)) + 4*vicon; % Q matrix is m11 in the document 
Q=double(Q);

% %% control _coefficient 
 c_num= 6*EIb*Eb*Ep*hb*d31*(hp+hb);
 c_den= (Eb^2*hb^4) + (4*Eb*Ep*hb^3*hp) + (6*Eb*Ep*hb^2*hp^2) + (4*Eb*Ep*hp^3*hb) + (Ep^2*hp^4); 
 c=c_num/c_den;
%  c= d31*EIb*(Eb*Ep*hb^2*hp + Eb*Ep*hb*hp^2)/(hp*((Eb^2*hb^4)/4 + (3*Ep*Eb*hb^3*hp)/4 + (Ep*Eb*hb^2*hp^2)/2));

% %% coefficients for tip deflection
% control_coefficient= vpa(c*(subs(du,lp_2)'- subs(du,lp_1)')); 
% tip_coefficients=subs(u,lb);
toc;

