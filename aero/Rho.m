function [atW,atT,atp,atrho,atVsonic]=Rho(z_in)   %z为海拔高度（m）
% 单位一律转为m
pSL=1.01325*10^5; %海平面大气压强
RhoSL=1.2250; %海平面大气密度
R0=6356.766;%地球半径
atW = zeros(size(z_in));
atT = zeros(size(z_in));
atp = zeros(size(z_in));
atrho = zeros(size(z_in));
atVsonic = zeros(size(z_in));
for i = 1:length(z_in)
    z = z_in(i)/1000;
    H=z/(1+z/R0);%地势高度
    if (z<0)
        T=0;
        p=0;
        Rho=0;
    end
    if(z>=0&&z<=11.0191)
        W=1-H/44.3308;
        T=288.15*W;
        p=pSL*W^5.2559;
        Rho=RhoSL*W^4.2559;
    end
    if(z>11.0191&&z<=20.0631)
        W=exp((14.9647-H)/6.3416);
        T=216.65;
        p=1.1953*10^-1*pSL*W;
        Rho=1.5898*10^-1*RhoSL*W;
    end
    if(z>20.0631&&z<=32.1619)
        W=1+(H-24.9021)/221.552;
        T=221.552*W;
        p=2.5158*10^-2*pSL*W^-34.1629;
        Rho=3.2722*10^-2*RhoSL*W^-25.1629;
    end
    if(z>32.1619&&z<=47.3501)
        W=1+(H-39.7499)/89.4107;
        T=250.350*W;
        p=2.8338*10^-3*pSL*W^-12.2011;
        Rho=3.2618*10^-3*RhoSL*W^-13.2011;
    end
    if(z>47.3501&&z<=51.4125)
        W=exp((48.6252-H)/7.9223);
        T=270.650;
        p=8.9155*10^-4*pSL*W;
        Rho=9.4920*10^-4*RhoSL*W;
    end
    if(z>51.4125&&z<=71.8020)
        W=1-(H-59.4390)/88.2218;
        T=247.021*W;
        p=2.1671*10^-4*pSL*W^12.2011;
        Rho=2.5280*10^-4*RhoSL*W^11.2011;
    end
    if(z>71.8020&&z<=86.0000)
        W=1-(H-78.0303)/100.2950;
        T=200.590*W;
        p=1.2274*10^-5*pSL*W^17.0816;
        Rho=1.7632*10^-5*RhoSL*W^16.0816;
    end
    if(z>86.0000)
        W=exp(87.2848-H)/5.4700;
        T=186.870*W;
        p=(2.2730+1.042*10^-3*H)*10^-6*W;
        Rho=3.6411*10^-6*W;
    end
    if(z<=91.0000)
        Vsonic=20.0468*sqrt(T);
    end
    if(z>91.0000)
        Vsonic=295;
    end
    
    atW(i) = W;
    atT(i) = T;
    atp(i) = p;
    atrho(i) = Rho;
    atVsonic(i) = Vsonic;
end
end
