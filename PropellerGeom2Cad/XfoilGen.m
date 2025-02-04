function [name] = XfoilGen(filePath,fileName,Thickness,TEgap,idx)
%XFOILGEN 이 함수의 요약 설명 위치
%   자세한 설명 위치
    clc
    mkdir 'Airfoil'\
    filePath=erase(filePath,' ');
    fileName=erase(fileName,' ');
    disp(" Airfoil Creating.......");
    tmpFile=split(fileName,'.');
    
    inpFileName=tmpFile{1}+ "_INP.inp";
    fid = fopen(inpFileName,'w');
    if (fid<=0)
        error([mfilename ':io'],'Unable to create xfoil.inp file');
        return;
    end
     fprintf(fid,'load %s\n',fileName);
    fprintf(fid,'\n\ngdes\n');
    %camber and Thickness
    fprintf(fid,'tset %f \n',Thickness);
    fprintf(fid,'\n');
     fprintf(fid,'eXec\ngset\n\n');


    fprintf(fid,'gdes\ntgap %f \n',TEgap/100.0);
    fprintf(fid,'\n');
    fprintf(fid,'eXec\ngset\n\n');

    fprintf(fid,'\nNAME %02d_%s_Based_%03d%03d\n',idx,tmpFile{1},floor(Thickness.*100),floor(TEgap.*100));
    fprintf(fid,'\nSAVE Airfoil\\%02d_%s_Based_%03d%03d.dat\n\n',idx,tmpFile{1},floor(Thickness.*100),floor(TEgap.*100));
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
    name=sprintf("%02d_%s_Based_%03d%03d.dat",idx,tmpFile{1},floor(Thickness.*100),floor(TEgap.*100));
end

