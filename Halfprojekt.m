clc
clear all
clf
format short

inrekrafter()
function D=inrekrafter()
Maxfart=90/3.6; %70 – 120 km/h
Maxacceleration=6; % 1a  6 m/s2
Maxretardation=-15; % 2a−  15 m/s2
Fordonetsvikt=180; %inkl. förare m  100 – 200 kg
g=9.82; %tyngdaccerelation
R=10; %Kurvradie 5-15m
c=0.25; %c  0,2–0,4 Luftmotstandskoefficient
A=0.6; % A  0,4 – 0,6 m2 Fordonetsfrontarea
df=950*10^-3; % Avstandfranframaxel till fordonets tyngdpunkt Fd  600 – 1000 mm
db=400*10^-3; % Avstandfranbakaxel till fordonets tyngdpunkt bd  200 – 400 mm
h=40*10^-2; % Tyngdpunktensvertikalapositionh 30 – 60 cm
h1=20*10^-2; %Avstandmellanluft tyngdpunkt och luft motstandets verdningslinje h1 15 – 25 cm
L=1150*10^-3; %Bakaxellangd L 1000 – 1200 mm
b1=100*10^-3; % Hjullagerposition1b  100 – 200 mm 
dh=300*10^-3; %Hjuldiameter 240 – 400 mm
rh=dh/2;
rb=dh*0.1; %Bromsbackens position 0.35 dh
rd=dh*0.1; %Drevets radie 0.4
bb=175*10^-3;
rho=1.21;
Karlradie=2;


Lf1=raktfram(Maxfart,Maxacceleration,Maxretardation,Fordonetsvikt,g,R,c,A,df,db,h,h1,L,b1,dh,rb,rd,bb,rho);
Lf2=acceleration(Maxfart,Maxacceleration,Maxretardation,Fordonetsvikt,g,R,c,A,df,db,h,h1,L,b1,dh,rb,rd,bb,rho);
Lf3=retardation(Maxfart,Maxacceleration,Maxretardation,Fordonetsvikt,g,R,c,A,df,db,h,h1,L,b1,dh,rb,rd,bb,rho);
Lf4=kurva(Maxfart,Maxacceleration,Maxretardation,Fordonetsvikt,g,R,c,A,df,db,h,h1,L,b1,dh,rb,rd,bb,rho);
Lf=[Lf1;Lf2;Lf3;Lf4];
%Snittningar
for i=1:4
inre=snittot(rh,rb,rd,L,b1,bb,Lf(i,1),Lf(i,2),Lf(i,3),Lf(i,4),Lf(i,5),Lf(i,6),Lf(i,7),Lf(i,8),Lf(i,9),Lf(i,10),Lf(i,11),Lf(i,12));
Mboj(i,:)=inre(1,:);
Tvar(i,:)=inre(2,:);
Mv(i,:)=inre(3,:);
Nx(i,:)=inre(4,:);
end


n=0.0000001;
x=(0:n:L);
figure(1)
plot(x,Mboj(:,:));hold on
title('$$\it M_{b}=\sqrt{\it M_{y}^{2}+M_{z}^{2}}$$',Interpreter='latex')
xlabel('$$\it X [m]$$',Interpreter='latex');ylabel('$$\it M_{b} [Nm]$$',Interpreter='latex');
legend('Maxfart', 'Maxacceleration', 'Maxretardation', 'Kurva');
hold off
figure(2)
plot(x,Tvar(:,:));hold on
title('$$ \it T=\sqrt{\it T_{y}^{2}+T_{z}^{2}}$$',Interpreter='latex')
xlabel('$$\it X [m]$$',Interpreter='latex');ylabel('$$\it T [N]$$',Interpreter='latex');
legend('Maxfart', 'Maxacceleration', 'Maxretardation', 'Kurva');
hold off
figure(3)
plot(x,Mv(:,:));hold on
title('$$\it M_{x}$$',Interpreter='latex')
xlabel('$$\it X [m]$$',Interpreter='latex');ylabel('$$\it M_{x}[Nm]$$',Interpreter='latex');
legend('Maxfart', 'Maxacceleration', 'Maxretardation', 'Kurva');
hold off
figure(4)
plot(x,Nx(:,:));hold on
title('$$\it N_{x}$$',Interpreter='latex')
xlabel('$$\it X [m]$$',Interpreter='latex');ylabel('$$\it N_{x} [N]$$',Interpreter='latex');
legend('Maxfart', 'Maxacceleration', 'Maxretardation', 'Kurva');
hold off
B1=(b1/L)*length(Mboj);
SFp=2.8;
Rp02=310*10^-3;


%För stora diametern 
MV=Mv(1:4,B1:end-B1);
Mb=Mboj(1:4,B1:end-B1);
NX=Nx(1:4,B1:end-B1);


Wv=@(D) (D^3)*pi/16; 
Wb=@(D) (D^3)*pi/32; 
Tvarsnitt=@(D) ((D/2)^2)*pi;
tau=@(D) (MV./Wv(D));
sigmatot= @(D) (NX./Tvarsnitt(D))+((Mb./Wb(D)));
effektivspanningar=@(D) (3*(tau(D).^2)+sigmatot(D).^2).^(1/2);
Maxeffektivspanning=@(D) (max(effektivspanningar(D), [], 'ALL'));
dim=@(D)((Rp02)/Maxeffektivspanning(D))-SFp;
D1=fzero(dim,50);


%För lilla diametern 
MV=[Mv(1:4,1:B1),Mv(1:4,end-B1:end)];
Mb=[Mboj(1:4,1:B1),Mboj(1:4,end-B1:end)];
NX=[Nx(1:4,1:B1),Nx(1:4,end-B1:end)];
Wv=@(D) (D^3)*pi/16; 
Wb=@(D) (D^3)*pi/32; 
Tvarsnitt=@(D) ((D/2)^2)*pi;
tau=@(D) (MV./Wv(D));
sigmatot= @(D) (NX./Tvarsnitt(D))+((Mb./Wb(D)));
effektivspanningar=@(D) (3*(tau(D).^2)+sigmatot(D).^2).^(1/2);
Maxeffektivspanning=@(D) (max(effektivspanningar(D), [], 'ALL'));
dim=@(D)((Rp02)/Maxeffektivspanning(D))-SFp;
D06=fzero(dim,50);


if D1<D06/0.6
    D=D06/0.6;
else D=D1;
end


%Normalspännig
%För stora diameter
MV=Mv;
Mb=Mboj;
NX=Nx;
Rp02=310;
Wv=@(D) (D^3)*pi/16; 
Wb=@(D) (D^3)*pi/32; 
Tvarsnitt=@(D) ((D/2)^2)*pi;
tau=@(D) (MV./Wv(D));
sigmatot= @(D) (NX./Tvarsnitt(D))+((Mb./Wb(D)));
effektivspanningar=@(D) (3*(tau(D).^2)+sigmatot(D).^2).^(1/2);
spannD=effektivspanningar(D);
sigmaD=sigmatot(D);
%För lilla diametern 
Tapp=D*0.6;
spann06=effektivspanningar(Tapp);
sigma06=sigmatot(Tapp);
spannplot=[spann06(1:4,1:B1),spannD(1:4,B1:end-B1-1),spann06(1:4,end-B1:end)];
sigmalplott=[sigma06(1:4,1:B1),sigmaD(1:4,B1:end-B1-1),sigma06(1:4,end-B1:end)];


% Spänningskoncentrationer vid lagren
qd=Karlradie/(D*0.6);
Ktv=1.6;
Ktb=2.0;
Ktn=2.3;
MV=[Mv(1:4,B1:B1),Mv(1:4,end-B1:end-B1)];
Mb=[Mboj(1:4,B1:B1),Mboj(1:4,end-B1:end-B1)];
NX=[Nx(1:4,B1:B1),Nx(1:4,end-B1:end-B1)];

Wv=@(D) (D^3)*pi/16; 
Wb=@(D) (D^3)*pi/32; 
Tvarsnitt=@(D) ((D/2)^2)*pi;
tau=@(D) (MV./Wv(D));
sigmatot= @(D) (Ktn*NX./Tvarsnitt(D))+((Ktb*Mb./Wb(D)));
effektivspanningar=@(D) (3*((Ktv*tau(D)).^2)+(sigmatot(D)).^2).^(1/2);
spanningskoncentrationer=effektivspanningar(D*0.6);
size(spanningskoncentrationer)
plottxkoncent=[x(B1:B1),x(end-B1:end-B1)];


%Effektivspännig
figure(5)
plot(x, spannplot*10^3); hold on;
xlabel('$$\it X [m]$$',Interpreter='latex');ylabel('$$\sigma_{e} [MPa]$$',Interpreter='latex');
plot(x,ones(1,length(x))*Rp02);
plot(plottxkoncent,spanningskoncentrationer*10^3,'*')
title('$$\sigma_{\it e}=\sqrt{\it\sigma_{\it e}^{2}+3\cdot\tau_{n}^{2}}$$',Interpreter='latex')
legend('Maxfart', 'Maxacceleration', 'Maxretardation', 'Kurva','sträckgräns','Maxfart koncentration', 'Maxacceleration koncentration', 'Maxretardation koncentration', 'Kurva koncentration');
hold off
%Normalspänning
figure(6)
plot(x,sigmalplott*10^3); hold on;
title('$$\sigma_{\it n}$$',Interpreter='latex')
xlabel('$$\it X [m]$$',Interpreter='latex'); ylabel('$$\sigma_{n} [MPa]$$',Interpreter='latex')
legend('Maxfart', 'Maxacceleration', 'Maxretardation', 'Kurva');
hold off
%Till highdiagram
maxismigma=max(sigmalplott, [], 'ALL')*10^3;

sigmaUb = 270; 
sigmaUbp = 240; 
Rm = 640;  
d = 0.6*D; 

lambda = 1; 
Kd = 1; 
Kr = 0.8;   
q = 0.8; 
Kt = (Karlradie/d)+(D/d); 
Kf = 1 + q*(Kt-1); 
nu = 1.9; 

Rf = lambda/(Kf*Kd*Kr); % Reduceringsfaktorn  

x = [0 sigmaUbp Rm];
y = [sigmaUb sigmaUbp 0];

% reducerade haj 

sigmaUb_red = sigmaUb*Rf; 
sigmaUbp_red = sigmaUbp*Rf;
x_red = [0 sigmaUbp Rm];
y_red = [sigmaUb_red sigmaUbp_red 0];
x_red2 = x_red; 
y_red2 = y_red./nu; 
figure(7)
plot(x,y)
hold on 
plot(x_red, y_red)
plot(x_red2, y_red2)
title('$$\it Haighdiagram$$',Interpreter='latex')
legend('Idial utmattning', 'Reduktion', 'Reduktion med säkerthetsfaktor')
text(0,maxismigma,'\bullet \leftarrow \sigma_{\itn}^{\itmax}');
xlabel('$$\it\sigma_{m} [MPa]$$', Interpreter='latex'); ylabel('$$\it\sigma_{a} [MPa]$$', Interpreter='latex')
xlim([0 650]); 
hold off
end
%_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
function Inre=snittot(rh,rb,rd,L,b1,bb,R1x,R1y,R1z,R2y,R2z,Fd,Vbi,Vby,Fk,Fb, Hbi,Hby)
n=0.0000001; 
x1=(0:n:b1); x2=x1(end):n:(((L/2)-bb));  x3=x2(end):n:((L/2)); 
x4=x3(end):n:(((L/2)+bb)); x5=x4(end):n:((L-b1)); x6=x5(end):n:(L);

inrekraft1=snitt1(x1,bb,Fd,rh,rb,Fb,Fk,rd,R1y,R2y,L,b1,Vbi,Hbi,R1z,R2z);
inrekraft2=snitt2(x2,bb,Fd,rh,rb,Fb,Fk,rd,R1y,R2y,L,b1,Vbi,Hbi,R1z,R2z);
inrekraft3=snitt3(x3,bb,Fd,rh,rb,Fb,Fk,rd,R1y,R2y,L,b1,Vbi,Hbi,R1z,R2z);
inrekraft4=snitt4(x4,bb,Fd,rh,rb,Fb,Fk,rd,R1y,R2y,L,b1,Vbi,Hbi,R1z,R2z);
inrekraft5=snitt5(x5,bb,Fd,rh,rb,Fb,Fk,rd,R1y,R2y,L,b1,Vbi,Hbi,R1z,R2z);
inrekraft6=snitt6(x6,bb,Fd,rh,rb,Fb,Fk,rd,R1y,R2y,L,b1,Vbi,Hbi,R1z,R2z);

IK=[inrekraft1,inrekraft2(:,2:end),inrekraft3(:,2:end),inrekraft4(:,2:end),inrekraft5(:,2:end),inrekraft6(:,2:end)];
formatspecT='%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s \n';
fprintf(formatspecT,'R1x','R1y','R1z','R2y','R2z','Fd','Vbi','Vby','Fk','Fb','Hbi','Hby')
formatspecn='%-15f %-15f %-15f %-15f %-15f %-15f %-15f %-15f %-15f %-15f %-15f %-15f \n';
fprintf(formatspecn,R1x,R1y,R1z,R2y,R2z,Fd,Vbi,Vby,Fk,Fb, Hbi,Hby)
Inre=[IK(1,:);IK(2,:);IK(3,:);IK(4,:)];

end
function inrekrafter=snitt1(x,~,Fd,rh,~,~,~,~,~,~,~,~,Vbi,Hbi,~,~)
% x=0 =>X=b1
Mh=rh*Fd/2;
My=-Vbi*x+Hbi*rh;
Mx=-Mh;
Mz=(-(Fd*x)/2);
Nx=Hbi;
Ty=-Fd/2;
Tz=-Vbi;
Mboj=sqrt((My.^2)+(Mz.^2));
Ttot=sqrt((Tz.^2)+(Ty.^2));
Mv=Mx;
inrekrafter=[Mboj;ones(1,length(x)).*Ttot;ones(1,length(x)).*Mv;ones(1,length(x)).*Nx];
end
function inrekrafter=snitt2(x,~,Fd,rh,~,~,~,~,R1y,~,~,b1,Vbi,Hbi,R1z,~)
% x=b1 =>X=(L/2)-bb
Mh=rh*Fd/2;
My=(-Vbi*x)-(R1z*(x-b1))+Hbi*rh;
Mx=-Mh;
Mz=(-Fd*x/2)-(R1y*(x-b1));
Nx=-Hbi;
Ty=-(Fd/2)-R1y;
Tz=-R1z-Vbi;
Mboj=sqrt((My.^2)+(Mz.^2));
Ttot=sqrt((Tz.^2)+(Ty.^2));
Mv=Mx;
inrekrafter=[Mboj;ones(1,length(x)).*Ttot;ones(1,length(x)).*Mv;ones(1,length(x)).*Nx];
end
function inrekrafter=snitt3(x,bb,Fd,rh,rb,Fb,~,~,R1y,~,L,b1,Vbi,Hbi,R1z,~)
% x=(L/2)-bb =>X=L/2
Mh=rh*Fd/2;
Mb=rb*Fb;
My=-Vbi*x-R1z*(x-b1)+Hbi*rh;
Mx=-Mh-Mb;
Mz=(-Fd*x/2)-R1y*(x-b1)+(Fb*(x-((L/2)-bb)));
Nx=-Hbi;
Ty=-(Fd/2)-R1y+Fb;
Tz=-R1z-Vbi;
Mboj=sqrt((My.^2)+(Mz.^2));
Ttot=sqrt((Tz.^2)+(Ty.^2));
Mv=Mx;
inrekrafter=[Mboj;ones(1,length(x)).*Ttot;ones(1,length(x)).*Mv;ones(1,length(x)).*Nx];
end
function inrekrafter=snitt4(x,bb,Fd,rh,rb,Fb,Fk,rd,R1y,~,L,b1,Vbi,Hbi,R1z,~)
% x=L/2 =>X=(L/2)+bb

Mh=rh*Fd/2;
Mb=rb*Fb;
Mk=Fk*rd;
My=-Vbi*x-R1z*(x-b1)+Hbi*rh;
Mx=-Mh-Mb+Mk;
Mz=(-Fd*x/2)-R1y*(x-b1)+(Fb*(x-((L/2)-bb)))-(Fk*(x-(L/2)));
Nx=-Hbi;
Ty=-(Fd/2)-R1y+Fb-Fk;
Tz=-R1z-Vbi;
Mboj=sqrt((My.^2)+(Mz.^2));
Ttot=sqrt((Tz.^2)+(Ty.^2));
Mv=Mx;
inrekrafter=[Mboj;ones(1,length(x)).*Ttot;ones(1,length(x)).*Mv;ones(1,length(x)).*Nx];
end
function inrekrafter=snitt5(x,bb,Fd,rh,rb,Fb,Fk,rd,R1y,~,L,b1,Vbi,Hbi,R1z,~)
% x=(L/2)+bb =>X=((L/2)-b1)
Mh=rh*Fd/2;
Mb=rb*Fb;
Mk=Fk*rd;
My=-Vbi*x-R1z*(x-b1)+Hbi*rh;
Mx=-Mh+Mk-2*Mb;
Mz=(-Fd*x/2)-R1y*(x-b1)+(Fb*(x-((L/2)-bb)))-(Fk*(x-(L/2)))+(Fb*(x-((L/2)+bb)));
Nx=-Hbi;
Ty=(-Fd/2)-R1y+2*Fb-Fk;
Tz=-R1z-Vbi;
Mboj=sqrt((My.^2)+(Mz.^2));
Ttot=sqrt((Tz.^2)+(Ty.^2));
Mv=Mx;
inrekrafter=[Mboj;ones(1,length(x)).*Ttot;ones(1,length(x)).*Mv;ones(1,length(x)).*Nx];
end
function inrekrafter=snitt6(x,bb,Fd,rh,rb,Fb,Fk,rd,R1y,R2y,L,b1,Vbi,Hbi,R1z,R2z)
%X=L-b1 => x=L
Mh=rh*Fd/2;
Mb=rb*Fb;
Mk=Fk*rd;
My=-Vbi*x-R1z*(x-b1)-R2z*(x-L+b1)+Hbi*rh;
Mx=-Mh+Mk-2*Mb;
Mz=(-Fd*x/2)-R1y*(x-b1)+(Fb*(x-((L/2)-bb)))-(Fk*(x-(L/2)))+(Fb*(x-((L/2)+bb)))-(R2y*(x-(L-b1)));
Nx=-Hbi;
Ty=-(Fd/2)-R1y +2*Fb-Fk-R2y;
Tz=-Vbi-R1z-R2z;
Mboj=sqrt((My.^2)+(Mz.^2));
Ttot=sqrt((Tz.^2)+(Ty.^2));
Mv=Mx;
inrekrafter=[Mboj;ones(1,length(x)).*Ttot;ones(1,length(x)).*Mv;ones(1,length(x)).*Nx];
end
%Lastfall
%_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
function yk=kurva(Maxfart,~,~,Fordonetsvikt,g,R,c,A,df,db,h,h1,L,b1,dh,~,rd,~,rho)
rh=dh/2;

Fl=(c*A*rho*(Maxfart.^2))/2;
Nb=((Fl*(h+h1))+(Fordonetsvikt*g*df))/(df+db);
Nf = (Fordonetsvikt*g)-Nb;
y=0.47;
fart=Maxfart*y;
Fc=Fordonetsvikt*((fart).^2)/R;
Fl=(c*A*rho*(fart.^2))/2;
Fd=Fl;
Fb=0;
a=0;
Fk=(Fd*rh)/(rd);
Vb = ((Fl*h1)+(Fordonetsvikt*g*df)+(Fd*h))/(df+db);
Vf = (Fordonetsvikt*g)-Vb;
Vbi = ((((Fordonetsvikt*g)/2)-((Fc*h)/L))*(Fl*(h+h1)+(Fordonetsvikt*g*df)))/(Fordonetsvikt*g*(df+db));
Vby = Nb-Vbi;

Hb = (Fc*df)/(db+df);
Hf = Fc - Hb;
Hbi = (Fc-Hf)/(1+(Vby/Vbi));
Hby = Hb - Hbi;

R1x = Hbi + Hby;
R2y = ((-Fd*(b1-(L/2)))-(Fk*(b1-(L/2))))/(L-(2*b1));
R1y = (Fd+Fk-R2y);
R1z = ((Vb*b1)-(Hb*dh*0.5)-(Vbi*L))./(L-(2*b1));
R2z = -R1z-Vb;
yk=[R1x,-R1y,R1z,-R2y,R2z,Fd,Vbi,Vby,Fk,Fb, Hbi,Hby];
end
function yk=retardation(Maxfart,~,Maxretardation,Fordonetsvikt,g,~,c,A,df,db,h,h1,L,b1,dh,rb,rd,bb,rho)
rh=dh/2;
Fl=(c*A*rho*(Maxfart.^2))/2;
Fd=Fl+(Maxretardation*Fordonetsvikt);
Fb=-(Fd*rh)/(2*rb);
Nb=((Fl*(h+h1))+(Fordonetsvikt*g*df)+(Fordonetsvikt*Maxretardation*h))/(df+db);
Fk=0;
R1x=0;
R2y=((Fb*2)-Fd)/2;
R1y=R2y;
R2z=-Nb/2;
R1z=R2z;
Vbi=Nb/2;
Vby=Nb/2;
Hbi=0;
Hby=0;
yk=[R1x,R1y,R1z,R2y,R2z,Fd,Vbi,Vby,Fk,Fb, Hbi,Hby];
end
function yk=acceleration(~,Maxacceleration,~,Fordonetsvikt,g,~,c,A,df,db,h,h1,L,b1,dh,~,rd,~,rho)
rh=dh/2;
Fd=(Fordonetsvikt*Maxacceleration);
Fb=0;
Nb=((Fordonetsvikt*g*df)+(Fordonetsvikt*Maxacceleration*h))/(df+db);
Fk=(Fd*rh)/(rd);
R1x=0;
R2y=(-Fd-Fk)/2;
R1y=R2y;
R2z=-Nb/2;
R1z=R2z;
Vbi=Nb/2;
Vby=Nb/2;
Hbi=0;
Hby=0;
yk=[R1x,R1y,R1z,R2y,R2z,Fd,Vbi,Vby,Fk,Fb, Hbi,Hby];
end
function yK=raktfram(Maxfart,~,~,Fordonetsvikt,g,~,c,A,df,db,h,h1,L,b1,dh,~,rd,~,rho)
rh=dh/2;
Fl=(c*A*rho*(Maxfart.^2))/2;
Fd=Fl;
Fb=0;
Nb=((Fl*(h+h1))+(Fordonetsvikt*g*df))/(df+db);
Fk=(Fd*rh)/(rd);
R1x=0;
R2y=(-Fd-Fk)/2;
R1y=R2y;
R2z=-Nb/2;
R1z=R2z;
Vbi=Nb/2;
Vby=Nb/2;
Hbi=0;
Hby=0;
yK=[R1x,R1y,R1z,R2y,R2z,Fd,Vbi,Vby,Fk,Fb, Hbi,Hby];
end








