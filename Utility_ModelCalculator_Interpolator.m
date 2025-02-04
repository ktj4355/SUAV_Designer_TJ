%Create Kriging model_Validation
%Cross Validation by K-fold Test
close all
clear
clc
CV_Cut=0.99;
%Section1 Model
CLtol=0.02;
CDtol=0.005;
iter=0;
listing=dir("Section1");
    listing(1);
    DB_Aero={}; %Alpha Cl Cd
    DBind=[]; %RE Mach
    Total_DB=[];
    for idx=3:size(listing,1)
        fileName=listing(idx).name;
        filePath=listing(idx).folder;
        fullName=[filePath '\' fileName];
        NameSplit=split(fileName,'_');
        AirfoilName=NameSplit{1};
        thick_major=str2num(NameSplit{2});
        thick_minor=str2num(NameSplit{3});
        thick=thick_major+thick_minor*0.01;
        

        %  alpha     CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge    XCp
        AeroDB_Local_line=readlines(fullName);
        Setting=double(split(AeroDB_Local_line(8)));
        Local_M=Setting(4);
        Local_Re=Setting(7)*(10^Setting(9));
        AeroDB_Local=double(split(AeroDB_Local_line(12:end-3)));
        AeroDB_Local(:,1)=[];
        DBind=[DBind;Local_Re,Local_M];
        DB_AeroLocal=[ones(size(AeroDB_Local,1),1)*thick ones(size(AeroDB_Local,1),1)*Local_Re  AeroDB_Local(:,1:3)];
        Total_DB=[Total_DB;DB_AeroLocal];
        DB_Aero=[DB_Aero;DB_AeroLocal];
    end

Section1_CL=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,4),"linear","nearest");
Section1_CD=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,5),"linear","nearest");
%% Section2
listing=dir("Section2");
    listing(1);
    DB_Aero={}; %Alpha Cl Cd
    DBind=[]; %RE Mach
    Total_DB=[];
    for idx=3:size(listing,1)
        fileName=listing(idx).name;
        filePath=listing(idx).folder;
        fullName=[filePath '\' fileName];
        NameSplit=split(fileName,'_');
        AirfoilName=NameSplit{1};
        thick_major=str2num(NameSplit{2});
        thick_minor=str2num(NameSplit{3});
        thick=thick_major+thick_minor*0.01;
        

        %  alpha     CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge    XCp
        AeroDB_Local_line=readlines(fullName);
        Setting=double(split(AeroDB_Local_line(8)));
        Local_M=Setting(4);
        Local_Re=Setting(7)*(10^Setting(9));
        AeroDB_Local=double(split(AeroDB_Local_line(12:end-3)));
        AeroDB_Local(:,1)=[];
        DBind=[DBind;Local_Re,Local_M];
        DB_AeroLocal=[ones(size(AeroDB_Local,1),1)*thick ones(size(AeroDB_Local,1),1)*Local_Re  AeroDB_Local(:,1:3)];
        Total_DB=[Total_DB;DB_AeroLocal];
        DB_Aero=[DB_Aero;DB_AeroLocal];
    end

Section2_CL=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,4),"linear","nearest");
Section2_CD=scatteredInterpolant(Total_DB(:,1),Total_DB(:,2),Total_DB(:,3),Total_DB(:,5),"linear","nearest");
ModelTitle=sprintf("InterpolatedModel");

save(ModelTitle,'Section1_CL','Section1_CD','Section2_CL','Section2_CD')



