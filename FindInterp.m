function Ri=FindInterp(Ti)
%This function is used for lookup and interpolation of given properties
%Created on 2020-4-1

%Original format Ri=FindInterp(melt.T,RLM,Ti)

global melt
global RLM

[dif,No]=min(abs(melt.T-Ti));%Ti is T interested

TC=[melt.T(No-1) melt.T(No) melt.T(No+1)];%nearest T
RC=[RLM(No-1) RLM(No) RLM(No+1)];%nearest property R
Ri=interp1(TC,RC,Ti,'spline');%cubic spline interpolation of property R

end

