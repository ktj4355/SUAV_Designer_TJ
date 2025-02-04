
%-------------------------------------------------------
%  SUAV Design Point Detecting Fucntion
%
%  Design Point
%
%------------------------------------------------------
clc
clear
close all
format longG
T=36; %N
D=2.5; %m
Vinf=12.69;
[Tmp, Pressure, den, D_vis, sonic] = STD_Atm(12000);
K_vis=D_vis./den;
Afan=((0.5*D).^2)*3.14;
n=480./60;
omega=n*2.*pi;

rR=0.1:0.05:1;
R=D./2;
B=2;
SpeedRatio=Vinf./(omega.*R);

%% 1) Select an initial Estimate for zeta
%Zeta : Displacement Velocity Ratio v'/V
ini_zeta = 0.01;

%% 2) Determine the values For F and phi  at Each Blade Section
zeta=ini_zeta;
zeta_old=100;
eta_old=100;
eta=0;
phi_tip=atan(SpeedRatio.*(1+zeta./2)); % radian
phi=atan(tan(phi_tip)./rR); %rad
phid=rad2deg(phi);

while(abs(eta-eta_old)>0.001)
    eta_old=eta;

    while(abs(zeta-zeta_old)>0.0001)
        eta_old=eta;

        phi_tip=phi(end); % radian
        phi_hub=phi(1); % radian



        f=(B./2).*(1-rR)./sin(phi_tip);
        Ftip=(2./pi).*acos(exp(-f));

        f=(B./2).*(rR-rR(1))./sin(phi_hub);
        Fhub=(2./pi).*acos(exp(-f));
        F=Ftip.*Fhub;
        F=Ftip; %Prantle Tip Correction
        %% 3) Determine the product Wc and Reynolds Number from EQ 16
        x=omega.*rR.*R.*Vinf;
        G=F.*x.*sin(phi).*cos(phi);
        Wc_Section=@(cl_section) 4.*pi.*SpeedRatio.*G.*Vinf.*R.*zeta./(cl_section.*B);
        cl_section=(0.8).*ones(1,size(rR,2));
        LD=(50).*ones(1,size(rR,2));
        e=1./LD;
        Wc=Wc_Section(cl_section) ;



        Re=den.*Wc./D_vis;




        %% 4) determind a and W
        
        a_ind=(zeta./2).*(cos(phi).^2).*(1-e.*tan(phi));
        a_swal=((zeta.*SpeedRatio)./(2.*rR)).*(cos(phi).*sin(phi)).*(1+e./(tan(phi)));
         a_ind_old=999;
         
        while(abs(a_ind-a_ind_old)>0.001)
            a_ind_old=a_ind;
            W_disk=omega.*(rR.*R).*(1-a_swal);
            
            W_n=Vinf.*(1+a_ind);
            W=sqrt(W_disk.^2+W_n.^2);
            c=Wc./W;
            Cy=cl_section.*(cos(phi)-e.*sin(phi));
            Cx=cl_section.*(sin(phi)+e.*cos(phi));

            Kind=Cy./(4.*sin(phi).*sin(phi));
            Kswrl=Cx./(4.*cos(phi).*sin(phi));

            sigma=B.*c./(2.*pi.*rR.*R);
            a_ind=sigma.*Kind./(F-sigma.*Kind);
            a_swal=sigma.*Kswrl./(F+sigma.*Kswrl);


        end

        %% 5) determind J and J
        dI1=4.*rR.*G.*(1-e.*tan(phi));
        dI2=SpeedRatio.*(dI1./(2.*rR)).*(1 + e./(tan(phi))).*sin(phi).*cos(phi);
        dJ1=4.*rR.*G.*(1+e./(tan(phi)));
        dJ2=            (dJ1./(2.*rR)).*(1 - e.*tan(phi))     .*cos(phi).*cos(phi);
        I1=0;
        I2=0;
        J1=0;
        J2=0;

        for i=1:length(rR)-1
            dr=(rR(i+1)-rR(i)).*R;
            I1=I1+0.5.*(dI1(i+1)+dI1(i)).*dr;
            I2=I2+0.5.*(dI2(i+1)+dI2(i)).*dr;

            J1=J1+0.5.*(dJ1(i+1)+dJ1(i)).*dr;
            J2=J2+0.5.*(dJ2(i+1)+dJ2(i)).*dr;
        end


        zeta_old=zeta;
        Tc_real=2.*T./(den.*(Vinf.^2).*pi.*(R.^2));
        zeta=(I1./(2.*I2))-sqrt(((I1./(2.*I2)).^2)-(Tc_real./I2));
        Pc=J1.*zeta+J2.*zeta.^2;
        Tc=I1.*zeta-I2.*zeta.^2;
        TCalc=(0.5*den.*(Vinf.^2).*pi.*(R.^2)).*Tc;
        PCalc=(0.5*den.*(Vinf.^3).*pi.*(R.^2)).*Pc;
        eta1=Tc./Pc;
        J=Vinf./(n.*D);
        QCalc=PCalc./(n*2*pi);
        Ct=TCalc./(den*(n^2)*(D^4));
        CP=PCalc./(den*(n^3)*(D^5));
        eta=(Ct./CP).*J;
    end
    %Estimation New Phi
    phi_tip=atan(SpeedRatio.*(1+zeta./2)); % radian
    phi=atan(Vinf.*(1+a_ind)./(omega.*(rR.*R).*(1-a_swal)));
    phi(end)=phi_tip;
    phid=rad2deg(phi);

end
eta; 
AR=R./mean(c);
figure(1); clf ;hold on
plot(rR, 0.5.*c./R)
plot(rR, -0.5.*c./R)
title("Designed SUAV Propeller Blade")
xlabel("r/R")
ylabel("c/R")
grid on
axis equal
figure(2); clf ;hold on
plot(rR, Re)
AOA=3.5;
beta= phid+AOA;
figure(3); clf ;hold on
plot(rR, phid)
plot(rR,beta)

TCalc=(0.5*den.*(Vinf.^2).*pi.*(R.^2)).*Tc;
PCalc=(0.5*den.*(Vinf.^3).*pi.*(R.^2)).*Pc;
Ct;
CP;
geom = [rR' (rR.*R)' c' beta'];

