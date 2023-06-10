function [nr,ni]=indwat(wl,xsal);
%
% from http://www.igf.fuw.edu.pl/meteo/stacja/kody/indwat.m
% Downloaded by Chuck DiMarzio, Northeastern University, July 2004
%  This does not handle arrays of wl correctly.
%
%  [nr,ni]=indwat(wl,xsal);
% Calculate water refration and extinction coefficient 
% input parameters:  wl=wavelength (in micrometers)
%                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
% output parameters: nr=index of refraction of sea water
%                    ni=extinction coefficient of sea water

% Indices of refraction for pure water from Hale and Querry, 
% Applied Optique, March 1973, Vol. 12,  No. 3, pp. 555-563
 twl=[0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,...
      0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,...
      0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,...
      1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,...
      2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,...
      3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,3.900,4.000];

 tnr=[1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,...
      1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,...
      1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,...
      1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,...
      1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,...
      1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,1.357,1.351];

 tni=[3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,...
      1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,...
      1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,...
      2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,...
      1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,...
      9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,...
      3.80E-03,4.60E-03];

 i=2;
 while i>0
    if (wl<twl(i))
      break
    end
    if i<62
     i=i+1;
    else
       break
    end   
  end            
            
 xwl=twl(i)-twl(i-1);        
 yr=tnr(i)-tnr(i-1);        
 yi=tni(i)-tni(i-1);        
 nr=tnr(i-1)+(wl-twl(i-1))*yr/xwl;
 ni=tni(i-1)+(wl-twl(i-1))*yi/xwl;
 
% Correction to be applied to the index of refraction and to the extinction 
% coefficients of the pure water to obtain the ocean water one (see for 
% example Friedman). By default, a typical sea water is assumed 
% (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
% In that case there is no correction for the extinction coefficient between 
% 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
% has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
% is a linear function of the salt concentration. Then, in 6S users are able 
% to enter the salt concentration (in ppt).
% REFERENCES:
% Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
% McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
%        New-York, 1965, p 129.
% Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
%        N.J., 1942, p 173.

  nrc=0.006;
  nic=0.000;
  nr=nr+nrc*(xsal/34.3);
  ni=ni+nic*(xsal/34.3);
  
