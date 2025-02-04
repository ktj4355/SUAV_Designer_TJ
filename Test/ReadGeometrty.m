clear 
close all
clc
% Read Geomety FIle
% r/R	r(m)	Chord (m)	Beta (deg)	AF_ThickRatio	AF

inputGeom=readmatrix("Geometry.xlsx");
geom=inputGeom(:,1:5);

%% SJDesign Calculator
alt=12000;
RPM=556;
Vinf=12.69;
B=2;
R=geom(2,end);
D=2*R;
SJ_data=Fucntion_BEMT_ROTATION(geom,alt,B,RPM,0,Vinf,0,0,0);
SJ_data_Result=SJ_data{5};
SJ_FullData=SJ_data{2};
SJ_R_Data=SJ_FullData{1};

figure(4); clf ;hold on

plot(SJ_R_Data(:,2)./R, SJ_R_Data(:,8))

%legend("Adkins", "BEMT")
clc
T_SJ_bemt=SJ_data_Result(2);
Q_SJ_bemt=SJ_data_Result(3);
P_SJ_bemt=SJ_data_Result(4);
n=RPM./60;
eta_SJ=(T_SJ_bemt.*Vinf)./(2.*pi.*n.*Q_SJ_bemt);