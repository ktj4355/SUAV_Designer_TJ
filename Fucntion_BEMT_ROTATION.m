function [calculatedData] = Fucntion_BEMT_ROTATION(PropGeom,alt,Blade,RPM,Vehicle_Speed_Foward,Vehicle_Speed_Asecnt,ForwardDirectionAngle,tilt_Angle,tilt_Phi)
%FUCNTION_BEMT_ROTATION 이 함수의 요약 설명 위치
%   자세한 설명 위치
% BEMT by Taejong Kim / APC tool
% Using APC Propeller Geomety
% Include Blended Section Analysis
% -> Require of 3 Section Interpolated Data

%Rotation Option Initial





%Load XLFR DAta
%uiopen("*.mat")
load("InterpolatedModel.mat");
%CL Model (RE, ALpha)


%Operating Condition
%[chord Beta] 1045 MR

R=PropGeom(end,2);
bmean=mean(PropGeom(:,3));
Rhub=0.8*0.0254;
AR=2*(R-Rhub)./bmean;


%r/R r Chord Beta(deg) TicknessRatio
BlendPosition=[2.85 4.65].*0.0254;
BlendPosition=BlendPosition./R;
BlendPosition=[1,1];
% Setup Condition

%alt     =0 ;         %   m
D	    =R.*2;       %   m
%V       =eps;
%   rev/min
%Blade   =2;

[Tmp, Pressure, rho, D_vis, a] = STD_Atm(alt);
%rho=0.3108

rho;
%Tmp=300;
%rho=1.18;
%rho=0.3108
n       =RPM./60;
%J	    =V./(n.*D);
V_tip	=2*pi*n.*R;
M_tip	=V_tip./340;
Afan    =pi*R*R;
V_Flow=sqrt(Vehicle_Speed_Foward^2+Vehicle_Speed_Asecnt^2);

OldT=0;
V_induced=0; % a of VLM
V_induced_old=0;
iter=0;
Iterdata=[];
Rotation_induced_factor=0; % a' of VLM
Rotation_induced_factor_old=0;




%% Calculate and Axis Convert Flow Vector

%Set Flow Vector
data_repeat=[];
calculatedData={};
idx=0;

idx=idx+1;
iter=0;

%Vehicle - Flow Axis Change
Vehicle_Direction_pol=[deg2rad(ForwardDirectionAngle),Vehicle_Speed_Foward,Vehicle_Speed_Asecnt];
diskNormal_shp=[deg2rad(tilt_Phi),deg2rad(90-tilt_Angle),1];
[V0x,V0y,V0z]=pol2cart(Vehicle_Direction_pol(1),Vehicle_Direction_pol(2),Vehicle_Direction_pol(3));
Vehicle_Direction=[V0x,V0y,V0z];
Inlet_Flow_Vector=-[V0x,V0y,V0z];

%Absolute Ground System Unit vector Setting
GroundUnitVector=[1 0 0;0 1 0; 0 0 1];                      %
GroundFlowVector=GroundUnitVector.*Inlet_Flow_Vector';      %

%Normal Vector Rotation Matrix R
pitch=deg2rad(tilt_Angle);
yaw=deg2rad(tilt_Phi);
Rp=[cos(pitch) 0 sin(pitch)
    0           1       0
    -sin(pitch) 0 cos(pitch)];
Ry=[cos(yaw) -sin(yaw) 0
    sin(yaw) cos(yaw) 0
    0 0 1];
Rot_mat=Ry*Rp;

% Disk Unit vector on Absolute Ground Frame
G_DiskUnitVector=Rot_mat*GroundUnitVector;
G_Flow_Axis_unit=G_DiskUnitVector(:,3);
G_Flow_TanX_unit=G_DiskUnitVector(:,1);
G_Flow_TanY_unit=G_DiskUnitVector(:,2);

% Disk Unit vector on Disk Frame
D_V_axis=dot(G_Flow_Axis_unit,Inlet_Flow_Vector);
D_V_tanX=dot(G_Flow_TanX_unit,Inlet_Flow_Vector);
D_V_tanY=dot(G_Flow_TanY_unit,Inlet_Flow_Vector);
D_FlowVector=[D_V_tanX D_V_tanY D_V_axis]';

% Disk Flow vector on Ground Frame
G_Flow_Axis_vec=D_V_axis.*G_Flow_Axis_unit;
G_Flow_TanX_vec=D_V_tanX.*G_Flow_TanX_unit;
G_Flow_TanY_vec=D_V_tanY.*G_Flow_TanY_unit;

%BLP : Blade Local Plane
FlowVector=D_FlowVector;

nAngle=36;
dAngle=2.*pi./nAngle;
BLP_Angle_SET=linspace(dAngle,2.*pi,nAngle);
J=abs(V0z)./(n.*D);
FlowData=[Vehicle_Speed_Foward Vehicle_Speed_Asecnt ForwardDirectionAngle tilt_Angle tilt_Phi];
% FlowData Vector Output Format
% 1) Vehicle_Speed_Foward         % 2) Vehicle_Speed_Asecnt
% 3) ForwardDirectionAngle        % 4) tilt_Angle
% 5) tilt_Phi

while(iter<100)
    data=[];
    T=0;
    Q=0;
    Traw=0;
    Qraw=0;

    iter = iter + 1;
    TotalForceVector=[0 0 0];
    %{
        
    figure(100); clf; hold on; grid on;
    txt=sprintf("Spanwise Blade AoA Distribution\nTilt Angle : %.1f deg | Forward Speed :%.1f m/s",tilt_Angle,Vehicle_Speed_Foward)
    xlabel("r/R");    ylabel("Blade AoA (deg)");    title(txt);

    figure(101); clf; hold on; grid on;
    txt=sprintf("Spanwise Blade Local Lift Coeefficient \nTilt Angle : %.1f deg | Forward Speed :%.1f m/s",tilt_Angle,Vehicle_Speed_Foward)
    xlabel("r/R");    ylabel("Lift Coefficient");    title(txt);
    figure(102); clf; hold on; grid on;
    txt=sprintf("Spanwise Blade Local Drag Coeefficient \nTilt Angle : %.1f deg | Forward Speed :%.1f m/s",tilt_Angle,Vehicle_Speed_Foward)
    xlabel("r/R");    ylabel("Drag Coefficient");    title(txt);
    %}
    iii=0;
    angleData=[];
    angleFlowData=[];
    angleData2=[];
    dataFull={};
    for angleindex=1:length(BLP_Angle_SET)
        tmpdata=[];
        data=[];



        % 1-1) Calculate Flow Vector
        BLP_Angle = BLP_Angle_SET(angleindex);
        BLP_NorVec_pol=[1,BLP_Angle,0]; %r the z
        BLP_NorVec_Car=[];%x y z
        [BLP_NorVec_Car(1),BLP_NorVec_Car(2),BLP_NorVec_Car(3)]=pol2cart(BLP_NorVec_pol(2),BLP_NorVec_pol(1),BLP_NorVec_pol(3)); %the r z
        BLP_TanVec_Car=[-BLP_NorVec_Car(2),BLP_NorVec_Car(1), 0];

        %
        BLP_inflow_Ax=dot(FlowVector,[0,0,-1]);
        BLP_inflow_Span=dot(FlowVector,BLP_NorVec_Car);
        BLP_inflow_Chord=dot(FlowVector,BLP_TanVec_Car);

        Vec_inflow_Ax=BLP_inflow_Ax.*[0,0,-1];
        Vec_inflow_Span=BLP_inflow_Span.*BLP_NorVec_Car;
        Vec_inflow_Chord=BLP_inflow_Chord.*BLP_TanVec_Car;

        BLP_Position=BLP_NorVec_Car.*R;
        angleFlowData=[angleFlowData ; BLP_Angle BLP_Position Vec_inflow_Ax Vec_inflow_Span Vec_inflow_Chord];
        LocalForceVector=[0 0 0];
        TdAngle=0;
        QdAngle=0;
        for r_idx=1:size(PropGeom,1)
            % Calculate Radius Section


            rR=PropGeom(r_idx,1);
            r=PropGeom(r_idx,2);
            if r_idx==1; dr=r;
            else;  dr=r-PropGeom(r_idx-1,2);
            end
            % Local Flow Calculate

            Vt0=2*pi*n*(1-Rotation_induced_factor).*r;   % Local Tanjential Velocity
            Vax=BLP_inflow_Ax+V_induced;
            Vt=Vt0-BLP_inflow_Chord;        % Local Axial Velocity
            Phi=rad2deg(atan2(Vax,Vt));
            % Local Flow angle
            vinVec=[BLP_inflow_Chord,-Vax]; % (Chord,Z)
            rotVec=[Vt0,0];
            displacementV=vinVec-rotVec;
            Vlocal=sqrt(Vax^2+Vt^2);
            b=PropGeom(r_idx,3);                        % Prop Chord
            Re=Vlocal.*rho.*b./D_vis;                   % Local Reynolds Number
            M_local=Vlocal./a;                          % Local Mach Number

            % Geometry Based AOA / CL Calculator
            beta=PropGeom(r_idx,4);                % Beta Geomety
            alpha=beta-Phi;                     % Local alpha
            ThickRatio=PropGeom(r_idx,5).*100;          % Geometry Thick Ratio
            % -> Calculating Real
            % Airfoil and Cl Cd
            % Claculate Cl Cd by Interpolation model
            inp_Cl=Section1_CL(ThickRatio, Re ,alpha);
            inp_Cd=Section1_CD(ThickRatio, Re ,alpha);

            %{
                if BlendPosition(1)>rR
                    airfoil=1;
                    inp_Cl=Section1_CL(ThickRatio, Re ,alpha);
                    inp_Cd=Section1_CD(ThickRatio, Re ,alpha);
                elseif BlendPosition(2)<rR
                    airfoil=2;
                    inp_Cl=Section2_CL(ThickRatio ,Re, alpha);
                    inp_Cd=Section2_CD(ThickRatio, Re, alpha);
                else
                    airfoil=0;
                    BR=abs((rR-BlendPosition(1))./(BlendPosition(2)-BlendPosition(1)));
                    inp_Cl=BLD_CL(BR, Re ,alpha);
                    inp_Cd=BLD_CD(BR ,Re, alpha);
                end

            %}

            CL0=Section1_CL(ThickRatio, Re ,0);
            inp_Cl=Section1_CL(ThickRatio, Re ,alpha);
            inp_Cd=Section1_CD(ThickRatio, Re ,alpha);
            Cl=inp_Cl;
            espan=1;
            inducedCD=((Cl-CL0).^2)./(pi.*espan.*AR);
            %inducedCD=0;
            CD=inp_Cd+inducedCD;
            %Calculate BET
            LD=Cl./CD;
            %Calculate BET
            gamma = rad2deg(atan(1/LD));
            K=Cl.*b./((sin(deg2rad(Phi)).^2).*cos(deg2rad(gamma)));
            Tc=K.*cos(deg2rad(Phi+gamma));
            Qc=r.*K.*sin(deg2rad(Phi+gamma));
            Fc=K.*sin(deg2rad(Phi+gamma));

            Phi_corr=Phi;
            if Phi<0;Phi_corr=eps;end
            % Prandtl Tip-Hub Correction


            Ptip=(Blade/2).*(R-r)./(r.*sin(deg2rad(Phi_corr)));
            Ftip=(2/pi)*acos(exp(-Ptip));
          %  Ftip=1;
            Phub=(Blade/2).*(r-Rhub)./(r.*sin(deg2rad(Phi_corr)));
            Fhub=(2/pi)*acos(exp(-Phub));
            Fhub=1;
            Tdrdth=(0.5*rho*Vax.^2).*Tc.*dr.*(1/nAngle);
            Qdrdth=(0.5*rho*Vax.^2).*Qc.*dr.*(1/nAngle);
            Fdrdth=(0.5*rho*Vax.^2).*Fc.*dr.*(1/nAngle);

            TdrCorr=Tdrdth*Ftip*Fhub;
            QdrCorr=Qdrdth*Ftip*Fhub;
            FdrCorr=Fdrdth*Ftip*Fhub;
            Tdr=[0,0,1].*TdrCorr;
            Fdr=-BLP_TanVec_Car.*FdrCorr;

            TdAngle=TdAngle+TdrCorr;
            QdAngle=QdAngle+QdrCorr;
            TotalForceVector=TotalForceVector+Tdr+Fdr;
            LocalForceVector=LocalForceVector+Tdr+Fdr;

            % data Vector Output Format
            % 1) Angle              % 2) r
            % 3) dr                 % 4) Vt
            % 5) Vax                % 6) Vlocal
            % 7) Phi                % 8) Cl
            % 9) Cd                 % 10) LD
            % 11) b (chord)         % 12) beta
            % 13) alpha             % 14) gamma
            % 15) Re                % 16) Mach
            % 17) dT (Correction)   % 18) dQ (Correction)   % 19) dF (Correction)
            % 20) dT (No Correction)% 21) dQ (No Correction)% 22) dF (No Correction)

            data=[data; BLP_Angle r dr Vt Vax Vlocal Phi Cl CD LD b beta alpha gamma Re M_local TdrCorr QdrCorr FdrCorr Tdrdth Qdrdth Fdrdth ] ;

        end
        dataFull(angleindex,1)={data};

        LocalForceVector=LocalForceVector *2;
        data_r=data(:,1)./R;
        data_a=data(:,11);
        data_cl=data(:,7);

        data_cd=data(:,21);
        T=T+TdAngle;
        Q=Q+QdAngle;
        iii=iii+1;
        %figure(111); yyaxis left     ;   hold on;        plot(data_r,data_a)
        %figure(100);  hold on;        plot(data_r,data_a,'k.-');
        %figure(101);  hold on;        plot(data_r,data_cl,'k.-');
        % figure(102);  hold on;        plot(data_r,data_cd,'k.-');
        %figure(111); yyaxis right     ;  hold on;        plot(data_r,data_cd)


        angleData=[angleData; BLP_Angle,LocalForceVector,QdAngle,TdAngle];
        % angle data Vector Output Format
        % 1) Angle              % 2) Local_dFx
        % 3) Local_dFy          % 4) Local_dT(Disk)
        % 5) Angle dQ           % 6) Angle dT


    end
    TotalForceVector=TotalForceVector.*2;

    T=T*Blade;
    Q=Q*Blade;
    P=n*2*pi*Q;
    Ct=T./(rho*(n^2)*(D^4));
    CQ=Q./(rho*(n^2)*(D^5));
    CP=P./(rho*(n^3)*(D^5));
    eta=(Ct./CP).*J;

    %Calculate Induced Velocity
    Vax0=-D_FlowVector(3);

    A=1;
    B=2;
    C=-2.*T./(rho.*Afan.*(Vax0.^2));
    bb=(-B+sqrt(B^2-4.*A.*C))./(2.*A);
    aa=(1-0.6084).*bb;
    V_induced=Vax0.*aa;
    Rotation_induced_factor_old=Rotation_induced_factor;
    Rotation_induced_factor=(Q.*pi)./(2.*(Afan.^2).*rho.*(Vax0.*(1+aa).*n));

    %Calculate Swirl deltaW

    Iterdata=[Iterdata;iter V_induced+Vax T Q P eta Rotation_induced_factor];
    % Iter data Vector Output Format
    % 1) iter               % 2) Induced Velocity
    % 3) T                  % 4) Q
    % 5) P                  % 6) eta
    % 7) Rotaltion Induced Factor

    if abs((V_induced-V_induced_old)./V_induced_old)<0.001
        if abs((Rotation_induced_factor-Rotation_induced_factor_old)./V_induced_old)<0.001
            %fileName=sprintf("SpanwiseAOA_Tilt_%d_Forward_%d.png",tilt_Angle,Vehicle_Speed_Foward); saveas(figure(100),fileName);
            % fileName=sprintf("SpanwiseCL_Tilt_%d_Forward_%d.png",tilt_Angle,Vehicle_Speed_Foward); saveas(figure(101),fileName);
            % fileName=sprintf("SpanwiseCD_Tilt_%d_Forward_%d.png",tilt_Angle,Vehicle_Speed_Foward); saveas(figure(102),fileName);

            break;
        end

    end
    V_induced_old=V_induced;
    OldT=T;
end
Global_Force_Axis=G_DiskUnitVector(:,3).*TotalForceVector(3);
Global_Force_TanX=G_DiskUnitVector(:,1).*TotalForceVector(1);
Global_Force_TanY=G_DiskUnitVector(:,2).*TotalForceVector(2);
GlobalForce=Global_Force_Axis+Global_Force_TanX+Global_Force_TanY;


data_repeat=[tilt_Angle T Q P GlobalForce'];
% data_repeat Vector Output Format
% 1) tilt_Angle         % 2) T_disk
% 3) Q_disk             % 4) Power
% 5) Fx (Global)        % 6) Fy (Global)    % 7) Fz (Global)

calculatedData(idx,:)={FlowData dataFull angleFlowData angleData data_repeat Iterdata};

end

