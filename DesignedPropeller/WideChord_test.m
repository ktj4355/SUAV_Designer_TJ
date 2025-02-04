
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
T=49%N
D=2.5; %m
Vinf=25.89;
alt=21000;
JRatio=0.65	;
n=Vinf./(JRatio.*D);
RPM=n*60
[Tmp, Pressure, den, D_vis, sonic] = STD_Atm(alt);
q=0.5*den*(Vinf.^2);

K_vis=D_vis./den;
Afan=((0.5*D).^2)*3.14;
n=RPM./60;
omega=n*2.*pi;
AOA=3.8;
rR=0.05:0.05:1;
R=D./2;
B=2;
SpeedRatio=Vinf./(omega.*R);
AOAdb=linspace(AOA,AOA,length(rR));
AOAdb(find(rR==0.85):end)=linspace(AOAdb(find(rR==0.85)),AOAdb(end),length(AOAdb(find(rR==0.85):end)));
AOAdb(1:max(find(rR<=0.2)))=linspace(0,AOAdb(max(find(rR<=0.2))),length(AOAdb(1:max(find(rR<=0.2)))));

thk=[0.12
    0.11
    0.11
    0.087
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074
    0.074];


%% 1) Select an initial Estimate for zeta
%Zeta : Displacement Velocity Ratio v'/V
ini_zeta = 1;

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

        f=(B./2).*(rR-0)./sin(phi_hub);
        Fhub=(2./pi).*acos(exp(-f));
        Fk=Ftip.*Fhub;
        F=Ftip;

       % F(end)=2*F(end-1)-F(end-2);
        %F=Ftip.*(Fhub+Fhub(end));
        %% 3) Determine the product Wc and Reynolds Number from EQ 16
                load("InterpolatedModel.mat");
        cl_section=0.8.*ones(1,size(rR,2));
                 x=omega.*rR.*R.*Vinf;
        G=F.*x.*sin(phi).*cos(phi);
        %G(:)=1
        Wc_Section=@(cl_section) 4.*pi.*SpeedRatio.*G.*Vinf.*R.*zeta./(cl_section.*B);

      % %{ 
            for(i=1:length(cl_section))
            Relocal=100;
            Relocal_old=0;
            while_idx=0;
              AOA=AOAdb(i);
            ThickRatio=thk(i)*100;
            while(abs(Relocal_old-Relocal)>0.01)
                while_idx=while_idx+1;
                Relocal_old=Relocal;
                local_Cl=Section1_CL(ThickRatio ,Relocal, AOA);
                local_Cl0=Section1_CL(ThickRatio ,Relocal, 0);

                local_Cd=Section1_CD(ThickRatio ,Relocal, AOA);
                %inducedCD=((local_Cl-local_Cl0).^2)./(pi.*espan.*AR);
                Wc_local=4.*pi.*SpeedRatio.*G(i).*Vinf.*R.*zeta./(local_Cl.*B);

                WcLocal=Wc_Section(local_Cl);

                Relocal=den.*Wc_local./D_vis;
                if (while_idx>1000)
                    Relocal=0.5*(Relocal+Relocal_old);
                    break;
                end
            end


       
        Re(i)=Relocal;
        cl_section(i)=Section1_CL(ThickRatio ,Relocal, AOA);
        LD(i)=cl_section(i)./Section1_CD(ThickRatio ,Relocal, AOA);
        end

%}        
       % cl_section=0.6.*ones(1,size(rR,2));

        Wc=Wc_Section(cl_section) ;
        Re=den.*Wc./D_vis;
      %  LD=10.*ones(1,size(rR,2));
        e=1./LD;



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

            c(1:4)=linspace(c(4)*0.75,c(4),4)
            c(end-1)=(2*c(end-2)-c(end-3));
            c(end)=(2*c(end-1)-c(end-2));
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
        Tc_real=T./(0.5*den.*(Vinf.^2).*pi.*(R.^2));
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
%c(1)=c(2);
%c(2)=(2*c(3)-c(4));
%c(1)=(2*c(2)-c(3));
%c(end-1)=(2*c(end-2)-c(end-3));

%c(end)=(2*c(end-1)-c(end-2));
%c(1:max(find(rR<=0.2)))=linspace(0.1,c(max(find(rR<=0.2))),length(c(1:max(find(rR<=0.2)))));

P1_up=plot(rR, 0.2*c./R,'b-');
P1_down=plot(rR, -0.75*c./R,'b-');

axis equal
figure(2); clf ;hold on

plot(rR, F)


beta= phid+AOAdb;
%beta(2)=beta(3);

beta(1)=beta(2);

figure(3); clf ;hold on
plot(rR, phid)
plot(rR,beta)

TCalc=(0.5*den.*(Vinf.^2).*pi.*(R.^2)).*Tc;
PCalc=(0.5*den.*(Vinf.^3).*pi.*(R.^2)).*Pc;
Ct;
CP;
geom = [rR' (rR.*R)' c' beta'];

geom = [rR' (rR.*R)' c' beta'];

geom(:,5)=thk;

%% Offdesign Calculator

alt=21000;
[Tmp, Pressure, den, D_vis, sonic] = STD_Atm(alt);
Vinf=sqrt(2.*q./den);
%Vinf=25.89;
  
n=Vinf./(JRatio.*D);
RPM=n*60;
RPM_design=RPM
%RPM=560
data=Fucntion_BEMT_ROTATION(geom,alt,B,RPM,0,Vinf,0,0,0);

data_Result=data{5};
FullData=data{2};
R_Data=FullData{1};

figure(4); clf ;hold on
title("local CL")
plot(rR,cl_section)
plot(R_Data(:,2)./R, R_Data(:,8))

legend("Adkins", "BEMT")
clc
T_bemt=data_Result(2);
Q_bemt=data_Result(3);
P_bemt=data_Result(4);
n=RPM./60;
eta=(T_bemt.*Vinf)./(2.*pi.*n.*Q_bemt);




%% SJDesign Calculator
JRatio=0.543	;  
inputGeom=readmatrix("Geometry_SJ.xlsx");
geom_sj=inputGeom(:,1:5);
n=Vinf./(JRatio.*D);
RPM=n*60;
RPM_SJ=RPM
B=2;
R=geom(end,2);
D=2*R;
SJ_data=Fucntion_BEMT_ROTATION(geom_sj,alt,B,RPM,0,Vinf,0,0,0);
SJ_data_Result=SJ_data{5};
SJ_FullData=SJ_data{2};
SJ_R_Data=SJ_FullData{1};
%%
figure(5); clf ;hold on
title("local Re")
xlabel("r/R")
ylabel("Re")
plot(R_Data(:,2)./R, R_Data(:,15))
plot(SJ_R_Data(:,2)./R, SJ_R_Data(:,15))
legend("designed","Sejong V8(J=0.548)")
clc
T_SJ_bemt=SJ_data_Result(2);
Q_SJ_bemt=SJ_data_Result(3);
P_SJ_bemt=SJ_data_Result(4);
n=RPM./60;
eta_SJ=(T_SJ_bemt.*Vinf)./(2.*pi.*n.*Q_SJ_bemt);
clc
fprintf("[Designed]   T:%.2f Q:%.2f P:%.2f Eta:%.4f  RPM : %.1f\n",T_bemt,Q_bemt,P_bemt,eta,RPM_design)
fprintf("[Refference] T:%.2f Q:%.2f P:%.2f Eta:%.4f  RPM : %.1f\n\n",T_SJ_bemt,Q_SJ_bemt,P_SJ_bemt,eta_SJ,RPM_SJ)

figure(1); hold on
P2_up=plot(geom_sj(:,1), 0.25*geom_sj(:,3)./R,'r-');
P2_down=plot(geom_sj(:,1), -0.75*geom_sj(:,3)./R,'r-');
legend([P1_up, P2_up],{"designed","Sejong V8"})

                     