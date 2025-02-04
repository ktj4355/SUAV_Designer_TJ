%% Airfoil Optimizer
% Regression and Kriging Model
% kim tae jong  |   ktj4355@gmail.com   |   010 4355 1390
% Sejong University |  Propulsion Aerodynamic Lab.

clc;clear;close all;
%% Defalt data
disp("Airfoil Desiginer by LHS Sampling")
disp("If you have Error, Please Check the L-H Sampling Data (1st Mode)")
mkdir 'Airfoil'\

delete("Airfoil\*.*")
aoa=-2:0.5:10;
Re=136992;
Mach=0.02351;
nCrit=9;
AnalyisisData=[];
%% Create Case Airfoile File

data=importdata("sampledata.mat");
sampleData=data;
CaseInd=1;
%insert base Coord data
[fileFullname,Wd] = uigetfile({'*.dat';'*.*'});
filename=split(fileFullname,".");
fileFullname=erase(fileFullname,' ');
Wd=erase(Wd,' ');
for sampleind=1:size(sampleData,1)
    CaseInd=sampleind;
    clc

    disp("# "+CaseInd+" Airfoil Creating......."+CaseInd+"/"+size(sampleData,1));
    inpFileName=filename(1)+ "_INP.inp";
    fid = fopen(inpFileName,'w');
    modifingDATA=sampleData(sampleind,:);

    if (fid<=0)
        error([mfilename ':io'],'Unable to create xfoil.inp file');
        continue;
    end
    fprintf(fid,'load %s\n',fileFullname);
    fprintf(fid,'\nppar\n');
    fprintf(fid,'N\n200\n');
    fprintf(fid,'\n\ngdes\n');
    %camber and Thickness
    fprintf(fid,'tset %f %f\n',modifingDATA(3)/100.0,modifingDATA(1)/100.0);
    %camber and Thickreness xCordinate
    fprintf(fid,'high %f %f\n',modifingDATA(4)/100.0,modifingDATA(2)/100.0);
    % Leading Edge
    fprintf(fid,'lera %f\n\n',modifingDATA(5));
    fprintf(fid,'eXec\ngset\n\n');

    fprintf(fid,'\nNAME %s\n',"case"+CaseInd);
    fprintf(fid,'\nSAVE %s\n',"ModifiedAirfoil\case"+CaseInd+".dat");
    fprintf(fid,'quit\n');

    fclose(fid);
    % input Xfoil.exe
    wd = fileparts(which(mfilename)); % working directory, where xfoil.exe needs to be
    cmd = sprintf('cd %s && xfoil.exe <%s> xfoil.out',wd,inpFileName);
    [status,result] = system(cmd);
    if (status~=0)
        disp(result);
        error([mfilename ':system'],'Xfoil execution failed! %s',cmd);
    end
  

end
 mkdir 'XFLR5 DATA'\
delete("XFLR5 DATA\*.*")

 input("Finish!, Press Enter to Mainmenu")
