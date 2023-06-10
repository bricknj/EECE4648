function [hbo2,hb]=hbspec(w);
% [k_hbo2,k_hb]=hbspec(w);
% Find Hb and HbO2 spectra for given wavelengths
%   w is a single number or a vector
% Wavelengths in list must be multiples of 2
%  in the range 250 to 1000
% T = 10^{\kappa C L}
data=load('2231.dat');
inds=find(ismember(data(:,1),w))';
hbo2=data(inds,2)';
hb=data(inds,3)';
if(length(w)~=length(hb));
   display('Warning: Use multiples of 2nm for wavelength, from 250 to 1000')
end;
