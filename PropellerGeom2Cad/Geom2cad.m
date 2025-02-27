%% Xfoil - Propeller Geometry Maker 
%% Taejong Kim
clc
clear
close all

%% Geom Data Read
% 같은 폴더안에 있어야함. Xlsx
inputGeom=readmatrix("Geom Test.xlsx");
geom=inputGeom(:,1:5);
if (geom(1,1)~=0)
geom=[geom(1,:);geom];
geom(1,1)=0;
geom(1,2)=0;
geom(1,3)=2*geom(2,3)-geom(3,3);
end
%% TEGAP Function
TEGAP=0.001*1./(geom(:,3));
TEGAP=TEGAP*100;
%% Airfoil Section Definitnion
mkdir 'Airfoil'\
delete("Airfoil\*.*")
% [Thickness TEgap]
Airfoil_Shape=[geom(:,5), TEGAP];
%% SectionPointAirfoil Create
[AF_filename, AF_file_path]=uigetfile(".dat");
name=[];
for(idx=1:size(Airfoil_Shape,1))
    nameTmp=XfoilGen(AF_file_path,AF_filename,Airfoil_Shape(idx,1),Airfoil_Shape(idx,2),idx);
    name=[name;nameTmp];
end

%% Rotate 0.25 Chord Point
clc
clf
airfoil_point={};
TE_line=[];
TE_line2=[];

LE_line=[];

for(idx=1:size(name,1))
    %for(idx=20)
    airfoilFilename=name(idx);
    A=readmatrix("Airfoil\"+airfoilFilename);

    sf=1000; %mm Scale Factor
    Scale_Chord=geom(idx,3)*sf;
    Scale_Span=geom(idx,2)*sf;
    T_scale=Scale_Chord.*[1 0; 0 1];
    A(:,1)=A(:,1)-0.25;

    B=A*T_scale;
    aoa=geom(idx,4);
    %aoa=10
    T_rotate=[cosd(aoa),sind(aoa);-sind(aoa) cosd(aoa)]; %반시계방향 회전 매트릭스
    C=(T_rotate*B')';

    figure(1)
    RealGeom=C;
    RealGeom(:,1)=Scale_Span; %x is Span
    RealGeom(:,2)=-C(:,1); %y is chord
    RealGeom(:,3)=C(:,2); %% Z is Disk
    hold on
    %plot(A(:,1), A(:,2))
    %plot(B(:,1), B(:,2))
    plot3(RealGeom(:,1), RealGeom(:,2),RealGeom(:,3),'k--')
    airfoil_point{idx,1}=RealGeom;

    TE_line(idx,:)=[RealGeom(1,1), RealGeom(1,2),RealGeom(1,3)];
    TE_line2(idx,:)=[RealGeom(end,1), RealGeom(end,2),RealGeom(end,3)];

    LE_line(idx,:)=[RealGeom(81,1), RealGeom(81,2),RealGeom(81,3)];

    figure(2)
    hold on
    plot(RealGeom(:,2), RealGeom(:,3))
    axis equal

end
figure(1)
hold on
plot3(TE_line(:,1), TE_line(:,2),TE_line(:,3),'r-')
plot3(TE_line2(:,1), TE_line2(:,2),TE_line2(:,3),'r-')

plot3(LE_line(:,1), LE_line(:,2),LE_line(:,3),'r-')
view(-45,-45)
axis equal


%% Output CatiaMacro File
%prop
fid=fopen("PropDesinner.csv",'w');
fprintf(fid,"StartLoft\n");
for idx=1:size(airfoil_point,1)
    localPoint=airfoil_point{idx};

    fprintf(fid,"StartCurve\n");
    for idxJ=1:size(localPoint,1)
        x=localPoint(idxJ,1);
        y=localPoint(idxJ,2);
        z=localPoint(idxJ,3);
        fprintf(fid,"%f,%f,%f\n",x,y,z);
    end
    fprintf(fid,"EndCurve\n");

end
fprintf(fid,"EndLoft\n");
fprintf(fid,"End\n");
fclose(fid)

%% LE
%
fid=fopen("LE.csv",'w');
fprintf(fid,"StartLoft\n");


    fprintf(fid,"StartCurve\n");
    for idxJ=1:size(LE_line,1)
        x=LE_line(idxJ,1);
        y=LE_line(idxJ,2);
        z=LE_line(idxJ,3);
        fprintf(fid,"%f,%f,%f\n",x,y,z);
    end
    fprintf(fid,"EndCurve\n");


fprintf(fid,"EndLoft\n");
fprintf(fid,"End\n");
fclose(fid)

%prop
fid=fopen("TE1.csv",'w');
fprintf(fid,"StartLoft\n");


    fprintf(fid,"StartCurve\n");
    for idxJ=1:size(TE_line,1)
        x=TE_line(idxJ,1);
        y=TE_line(idxJ,2);
        z=TE_line(idxJ,3);
        fprintf(fid,"%f,%f,%f\n",x,y,z);
    end
    fprintf(fid,"EndCurve\n");


fprintf(fid,"EndLoft\n");
fprintf(fid,"End\n");
fclose(fid)


fid=fopen("TE2.csv",'w');
fprintf(fid,"StartLoft\n");


fprintf(fid,"StartCurve\n");
for idxJ=1:size(TE_line2,1)
    x=TE_line2(idxJ,1);
    y=TE_line2(idxJ,2);
    z=TE_line2(idxJ,3);
    fprintf(fid,"%f,%f,%f\n",x,y,z);
end
fprintf(fid,"EndCurve\n");


fprintf(fid,"EndLoft\n");
fprintf(fid,"End\n");
fclose(fid)
