
clc
close all
clear
disp("Interpolated Data...")
uiopen("*.mat")
%Sec1_KrigingData={dmodel_sec1 perf_sec1 CV_rmse S Y S_Validate Y_Validate}
%%
PropGeom=...
    [0.16	0.02032	0.01867408	34.2646	0.1557
    0.172	0.021844	0.01977644	34.8944	0.1443
    0.184	0.023368	0.02081276	34.9951	0.1342
    0.196	0.024892	0.0217805	34.6686	0.1253
    0.208	0.026416	0.0226822	33.9804	0.1176
    0.22	0.02794	0.02351532	32.9732	0.1109
    0.232	0.029464	0.02427986	31.6917	0.1052
    0.24502	0.03111754	0.02503678	30.3101	0.1
    0.2678	0.0340106	0.02616962	28.141	0.0931
    0.2916	0.0370332	0.02709926	26.1611	0.0882
    0.3154	0.0400558	0.02776982	24.4253	0.0851
    0.3392	0.0430784	0.02818384	22.8937	0.0831
    0.363	0.046101	0.02834386	21.5341	0.0816
    0.3868	0.0491236	0.02824988	20.3206	0.0801
    0.4106	0.0521462	0.02801112	19.2315	0.0788
    0.4344	0.0551688	0.0277114	18.2495	0.0776
    0.4582	0.0581914	0.0273558	17.3599	0.0765
    0.482	0.061214	0.02694432	16.5508	0.0756
    0.5058	0.0642366	0.0264795	15.8118	0.0748
    0.5296	0.0672592	0.02596388	15.1345	0.0741
    0.5534	0.0702818	0.02539746	14.5117	0.0735
    0.5772	0.0733044	0.02478024	13.9371	0.073
    0.601	0.076327	0.02411984	13.4055	0.0726
    0.6248	0.0793496	0.02341372	12.9123	0.0722
    0.6486	0.0823722	0.02266442	12.4535	0.072
    0.6724	0.0853948	0.02187702	12.0258	0.0718
    0.6962	0.0884174	0.02104898	11.6261	0.0717
    0.72	0.09144	0.02018284	11.2517	0.0717
    0.7438	0.0944626	0.01928368	10.9005	0.0717
    0.7676	0.0974852	0.01834896	10.5702	0.0718
    0.7914	0.1005078	0.0173863	10.2592	0.0719
    0.8152	0.1035304	0.01639062	9.9658	0.072
    0.839	0.106553	0.01536954	9.6885	0.0722
    0.8628	0.1095756	0.01432306	9.4261	0.0724
    0.8866	0.1125982	0.01325118	9.1774	0.0726
    0.9104	0.1156208	0.01215644	8.9414	0.0729
    0.9342	0.1186434	0.01104392	8.7172	0.0731
    0.95784	0.12164568	0.0097282	8.5052	0.0734
    0.98018	0.12448286	0.00712724	8.3142	0.0736
    1	0.127	0.00000254	8.1839	0.0738];

R=PropGeom(end,2);
BlendPosition=[2.85 4.65].*0.0254;
BlendPosition=BlendPosition./R;

data=[];
Xdata=[];
for ii=1:size(PropGeom,1)
    rR=PropGeom(ii,1);
    r=PropGeom(ii,2);
    chord=PropGeom(ii,3);
    beta=PropGeom(ii,4);
    ThickRatio=PropGeom(ii,5).*100;
    alpha=0;
    Re=75000;
    M=0.15;

    if BlendPosition(1)>rR
        airfoil=1;
        X=[ThickRatio Re alpha];
        inp_Cl=Section1_CL(ThickRatio, Re ,alpha);
        inp_Cd=Section1_CD(ThickRatio, Re ,alpha);
        out=[inp_Cl inp_Cd];
    elseif BlendPosition(2)<rR
        airfoil=2;
        X=[ThickRatio Re alpha];

        inp_Cl=Section2_CL(ThickRatio ,Re, alpha);
        inp_Cd=Section2_CD(ThickRatio, Re, alpha);
        out=[inp_Cl inp_Cd];
    else
        airfoil=0;
        BR=abs((rR-BlendPosition(1))./(BlendPosition(2)-BlendPosition(1)));
        X=[BR Re alpha];
   
        inp_Cl=BLD_CL(BR, Re ,alpha);
        inp_Cd=BLD_CD(BR ,Re, alpha);
        out=[inp_Cl inp_Cd];
    end

    cl0=out(1);
    cd0=out(2);
    ld0=cl0./cd0;
    gamma=rad2deg(atan(1/ld0));
    data=[data;rR airfoil cl0 cd0 ld0 gamma];
    Xdata=[Xdata;X];
end
plot(data(:,1),data(:,3))