clc
clear
close all
format longG
Vinf=36.19;
alt=25000;
Vinf=25.89;
alt=21000;
%Vinf=12.69;
%alt=12000;
D=2.5; %m
R=D./2;
B=2;
RPM_data=[]
J_Designed=[];
T_Designed=[];
Q_Designed=[];
eta_Designed=[];

J_Ref=[];
T_Ref=[];
Q_Ref=[];
eta_Ref=[];

J_Lib=0.1:0.05:0.7;
targetJ=0.5438
JTargetInd=max(find(J_Lib<0.537))
J_Lib=[J_Lib(1:JTargetInd) targetJ J_Lib(JTargetInd+1:end)];
targetJ=0.635
JTargetInd=max(find(J_Lib<0.635))
J_Lib=[J_Lib(1:JTargetInd) targetJ J_Lib(JTargetInd+1:end)];
for J_ind=1:length(J_Lib)

inputGeom=readmatrix("Geom_250110_J60.xlsx");
geom=inputGeom(:,1:5);
rR=geom(:,1);


[Tmp, Pressure, den, D_vis, sonic] = STD_Atm(alt);
q=0.5*den*(Vinf.^2);

%Vinf=sqrt(2.*q./den)
JRatio=J_Lib(J_ind);  


n=Vinf./(JRatio.*D);
RPM=n*60;
%RPM=560
data=Fucntion_BEMT_ROTATION(geom,alt,B,RPM,0,Vinf,0,0,0);

data_Result=data{5};
FullData=data{2};
R_Data=FullData{1};

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

T_SJ_bemt=SJ_data_Result(2);
Q_SJ_bemt=SJ_data_Result(3);
P_SJ_bemt=SJ_data_Result(4);
n=RPM./60;
eta_SJ=(T_SJ_bemt.*Vinf)./(2.*pi.*n.*Q_SJ_bemt);
fprintf("[J = %.2f (RPM = %d)] \n  ",JRatio,floor(RPM))
fprintf("[Designed]   T:%.2f Q:%.2f P:%.2f Eta:%.4f\n",T_bemt,Q_bemt,P_bemt,eta)
fprintf("[Refference] T:%.2f Q:%.2f P:%.2f Eta:%.4f\n\n",T_SJ_bemt,Q_SJ_bemt,P_SJ_bemt,eta_SJ)
J_Designed(J_ind)=JRatio;
T_Designed(J_ind)=T_bemt;
Q_Designed(J_ind)=Q_bemt;
eta_Designed(J_ind)=eta;

J_Ref(J_ind)=JRatio;
T_Ref(J_ind)=T_SJ_bemt;
Q_Ref(J_ind)=Q_SJ_bemt;
eta_Ref(J_ind)=eta_SJ;
RPM_data(J_ind)=RPM;
end

%%
figure(1)
clf
hold on 
plot(J_Designed,eta_Designed,'rx-');
plot(J_Ref,eta_Ref,'bx-');
legend("Designed","Sejong V8",Location='best')
title("Efficiency")
grid on

figure(2)
clf
hold on 
plot(J_Designed,Q_Designed,'rx-');
plot(J_Ref,Q_Ref,'bx-');
legend("Designed","Sejong V8",Location='best')
title("Torque")

grid on

figure(3)
clf
hold on 
plot(J_Designed,T_Designed,'rx-');
plot(J_Ref,T_Ref,'bx-');
legend("Designed","Sejong V8",Location='best')
title("Thrust")

grid on
figure(4)
clf
hold on 
plot(RPM_data,eta_Designed,'rx-');
plot(RPM_data,eta_Ref,'bx-');
legend("Designed","Sejong V8",Location='best')
title("Efficiency (RPM)")
grid on


figure(5)
clf
hold on 
plot(T_Designed,eta_Designed,'rx-');
plot(T_Ref,eta_Ref,'bx-');
title("Efficiency by Thrust")
xline(32)
legend("Designed","Sejong V8","32N line",Location='best')

xlabel("Thrust")
ylabel("Efficiency")
grid on

figure(6)
clf
hold on 
plot(T_Designed,Q_Designed,'rx-');
plot(T_Ref,Q_Ref,'bx-');
title("Efficiency by Thrust")
xline(32)
legend("Designed","Sejong V8","32N line",Location='best')

xlabel("Thrust")
ylabel("Torque")
grid on