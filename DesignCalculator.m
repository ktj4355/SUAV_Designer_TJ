clc
clear
close all
format longG
  Vinf=36.19;
  alt=25000;

  Vinf=25.89;
  alt=21000;

  Vinf=12.69;
  alt=12000;

D=2.5; %m
R=D./2;
B=2;
inputGeom=readmatrix("12km_HighEff.xlsx");
geom=inputGeom(:,1:5);
rR=geom(:,1);


[Tmp, Pressure, den, D_vis, sonic] = STD_Atm(alt);
q=0.5*den*(Vinf.^2);

%Vinf=sqrt(2.*q./den)
JRatio=0.543857142857143;  


n=Vinf./(JRatio.*D);
RPM=n*60;
%RPM=560
data=Fucntion_BEMT_ROTATION(geom,alt,B,RPM,0,Vinf,0,0,0);

data_Result=data{5};
FullData=data{2};
R_Data=FullData{1};

figure(4); clf ;hold on
title("local CL")
plot(R_Data(:,2)./R, R_Data(:,8))

legend("BEMT")
clc
T_bemt=data_Result(2);
Q_bemt=data_Result(3);
P_bemt=data_Result(4);
n=RPM./60;
eta=(T_bemt.*Vinf)./(2.*pi.*n.*Q_bemt);




%% SJDesign Calculator
inputGeom=readmatrix("Geometry_SJ.xlsx");
geom_sj=inputGeom(:,1:5);
%n=Vinf./(JRatio.*D);
%RPM=n*60;
%RPM=360

B=2;
R=geom_sj(end,2);
D=2*R;
SJ_data=Fucntion_BEMT_ROTATION(geom_sj,alt,B,RPM,0,Vinf,0,0,0);
SJ_data_Result=SJ_data{5};
SJ_FullData=SJ_data{2};
SJ_R_Data=SJ_FullData{1};
%%
figure(5); clf ;hold on;grid on

plot(R_Data(:,2)./R, R_Data(:,15),'r-')
%plot(SJ_R_Data(:,2)./R, SJ_R_Data(:,15),'b-')
legend("designed","Sejong V8");grid on
figure(6); clf ;hold on;grid on

plot(R_Data(:,2)./R, R_Data(:,13),'r-')
plot(SJ_R_Data(:,2)./R, SJ_R_Data(:,13),'b-')
legend("designed","Sejong V8");grid on
clc
T_SJ_bemt=SJ_data_Result(2);
Q_SJ_bemt=SJ_data_Result(3);
P_SJ_bemt=SJ_data_Result(4);
n=RPM./60;
eta_SJ=(T_SJ_bemt.*Vinf)./(2.*pi.*n.*Q_SJ_bemt);
clc
fprintf("[Designed]   T:%.4f Q:%.4f P:%.4f Eta:%.4f\n",T_bemt,Q_bemt,P_bemt,eta)
fprintf("[Refference] T:%.4f Q:%.4f P:%.4f Eta:%.4f\n",T_SJ_bemt,Q_SJ_bemt,P_SJ_bemt,eta_SJ)

figure(1); hold on
P1_up=plot(geom(:,1), 0.25*geom(:,3)./R,'r-');
P1_down=plot(geom(:,1), -0.75*geom(:,3)./R,'r-');
P2_up=plot(geom_sj(:,1), 0.25*geom_sj(:,3)./R,'b-');
P2_down=plot(geom_sj(:,1), -0.75*geom_sj(:,3)./R,'b-');
legend([P1_up, P2_up],{"designed","Sejong V8"})
axis equal
