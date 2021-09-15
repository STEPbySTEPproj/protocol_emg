function ic = cocontraction_winter( time, data1, data2)        
%COCONTRACTION Summary of this function goes here
%   Detailed explanation goes here


sumIntegrals = discrete_integrate(time, data1) + discrete_integrate(time, data2);

overlapCurve = min(data1,data2);
integralOverlap = discrete_integrate(time, overlapCurve);

ic = 2 * integralOverlap/sumIntegrals;


