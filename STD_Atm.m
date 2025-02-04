function [Tmp, Pressure, den, D_vis, a] = STD_Atm(Height)
%STD_ATM 이 함수의 요약 설명 위치
%   자세한 설명 위치
if Height <11000
    Tmp=288.15-0.0065.*Height;
    Pressure=101.325.*(Tmp./288.15).^5.2559;
elseif Height<25000
    Tmp=273.15-56.46;
    Pressure=22.65*exp(1.73-0.000157.*Height);
elseif Height>=25000
    Tmp=273.15-131.21+0.00299*Height;
    Pressure=2.488.*(Tmp./216.15).^(-11.388);
end

den=Pressure./(0.2869.*Tmp);
D_vis=((0.000001458).*(Tmp.^(1.5)))/(Tmp+110.4);
a=sqrt(1.4.*287.*Tmp);
end