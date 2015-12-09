function xtext=xtitles(iforw,iparm)

% This function finds the variable name as it has been used by saga
% PG, july 95
% xtext = xtitles(iforw,iparm)
% iforw = 1,2,7:  SNAP and SNAPRD and PROSIM
% iforw = 3,4  : OASES
% iforw = 5    : POPP
% iforw = 6    : TPEM
% iforw = 9    : RAMGEO
% iforw = 10   : GAMA
% iforw = 11   : ORCA
%
% To organize variable names that consist of different sizes,
% the cell array is used.  .cfh. 

if (iforw==1 | iforw==2 | iforw==7 | iforw==8)
%*** create text strings for physical parameters
%                   123456789012345678901234567890
% This part is for SNAP and SNAPRD and PROSIM
      phystxt{1} = 'Water depth (m)';
      phystxt{2} = 'Water sound speed (m/s)    ';
      phystxt{3} = 'Sediment sound speed (m/s) ';
      phystxt{4} = 'Attenuation (dB/\lambda)    ';
      phystxt{5} = 'Surface roughness (m)      ';
      phystxt{6} = 'Sediment density (g/cm^3)   ';
      phystxt{8} = 'Source depth (m)           ';
      phystxt{9} = 'Source range (km)          ';
      phystxt{11}= 'Shape coefficient          ';
      phystxt{12}= 'Bottom sound speed (m/s)   ';
      phystxt{13}= 'Bottom S-speed (m/s)       ';
      phystxt{14}= 'Bottom density (g/cm^3)     ';
      phystxt{15}= 'Receiver depth (m)         ';
      phystxt{16}= 'Depth of water speed pt.(m)';
      phystxt{17}= 'Sediment depth (m)         ';
      phystxt{18}= 'Depth of sed. speed pt.(m) ';
      phystxt{19}= 'Array tilt (m)             ';
      phystxt{20}= 'Sound speed profile #      ';
      phystxt{21}= 'Out of plane parameter     ';
      phystxt{22}= 'Model source streght (dB)  ';
      phystxt{23}= 'Second source param        ';
      phystxt{24}= 'Mean grain size (\phi)     ';
      phystxt{25}= 'Range offset (m)     ';
      phystxt{26}= 'Time offset (s)     ';
if (iforw==7)
      phystxt{20}= 'Time delay                 ';
      phystxt{21}= 'Horizontal tilt (deg)      ';

end
elseif (iforw==3 | iforw==4)
% This part is for OASES
      phystxt{1} = 'Depth of layer (m)         ';
      phystxt{2} = 'P-sound speed (m/s)        ';
      phystxt{3} = 'S-sound speed (m/s)        ';
      phystxt{4} = 'P-Attenuation (dB/\lambda)  ';
      phystxt{5} = 'S-Attenuation (dB/\lambda)  ';
      phystxt{6} = 'Sediment density (kg/m^3)   ';
      phystxt{7} = 'Layer thickness (m)        ';
      phystxt{8} = 'Source depth (m)           ';
      phystxt{9} = 'Source-receiver range (km) ';
      phystxt{11}= 'Shape coefficient          ';
      phystxt{15}= 'Receiver depth (m)         ';
elseif (iforw==5)
%*** create text strings for physical parameters
%                   123456789012345678901234567890
% This part is for a POPP
      phystxt{1} = 'Water depth (m)            ';
      phystxt{2} = 'Water sound speed (m/s)    ';
      phystxt{3} = 'Sediment sound speed (m/s) ';
      phystxt{4} = 'Attenuation sed (dB/\lambda)';
      phystxt{5} = 'sediment thickness (m)     ';
      phystxt{6} = 'Sediment density (g/cm^3)   ';
      phystxt{8} = 'Source depth (m)           ';
      phystxt{9} = 'Source range (km)          ';
      phystxt{11}= 'Shape coefficient          ';
      phystxt{15}= 'Receiver depth (m)         ';
      phystxt{16}= 'Depth of water speed pt.(m)';
      phystxt{17}= 'Lambert coefficient        ';
      phystxt{18}= 'Lambert power              ';
%***  
elseif (iforw==6)                                   % TPEM
      phystxt{1} = 'Refractivity               ';
      phystxt{2} = 'Heights of refract. pt (m) ';
      phystxt{3} = 'Terrain x-coordinates (m)  ';
      phystxt{4} = 'Terrain y-coordinates (m)  ';
      phystxt{5} = 'Source-receiver dist (m)  ';
      phystxt{6} = 'Source height (m)          ';
      phystxt{9} = 'Source-receiver dist  (km) ';
      phystxt{11}= 'Shape coefficient          ';
      phystxt{12}= 'Baseheight (m)             ';
      phystxt{13}= 'Thickness (m)              ';
      phystxt{14}= 'Offset (M-units)           ';
      phystxt{15}= 'M-defecit (M-units)        ';
      phystxt{16}= 'Noise-floor                ';
      phystxt{17}= 'Slopes ......              ';
      phystxt{18}= 'Clutter cross section (dB) ';
      phystxt{19}= 'Range factors   baseheight ';
  elseif (iforw==9)                                   % RAMGEO
      phystxt{1} = 'Ocean sound speed (m/s)   ';
      phystxt{2} = 'Bottom sound speed (m/s)  ';
      phystxt{3} = 'Bottom density (g/cm^3)   ';
      phystxt{4} = 'Bottom attenuation (dB)   ';
      phystxt{5} = 'Depth point (m)           ';
      phystxt{6} = 'Depth point (m)           ';
      phystxt{7} = 'Depth point (m)           ';
      phystxt{8} = 'Source Depth (m)          ';
      phystxt{9} = 'Source-receiver range (m) ';
      phystxt{11}= 'Shape coefficient         ';
      phystxt{10}= 'Depth point (m)           ';
      phystxt{12}= 'Bathymetry  (m)           ';
      phystxt{13}= 'Bathymetry range point (m)';
elseif (iforw==10 )
%GAMA
      phystxt{1} = 'Water depth (m)            ';
      phystxt{2} = 'Water sound speed (m/s)    ';
      phystxt{3} = 'Sediment sound speed (m/s) ';
      phystxt{4} = 'Attenuation (dB/\lambda)    ';
      phystxt{5} = '                           ';
      phystxt{6} = 'Sediment density (g/cm^3)   ';
      phystxt{8} = 'Source depth (m)           ';
      phystxt{9} = 'Source range (km)          ';
      phystxt{11}= 'Shape coefficient          ';
      phystxt{12}= 'Bottom sound speed (m/s)   ';
      phystxt{13}= 'Bottom S-speed (m/s)       ';
      phystxt{14}= 'Bottom density (g/cm^3)     ';
      phystxt{15}= 'Receiver depth (m)         ';
      phystxt{16}= 'Depth of water speed pt.(m)';
      phystxt{17}= 'Sediment depth (m)         ';
      phystxt{18}= 'Depth of sed. speed pt.(m) ';
      phystxt{19}= 'Array tilt (m)             ';
      phystxt{20}= 'Sound speed profile #      ';
      phystxt{21}= ' Horizontal Array tilt (m) ';
      phystxt{22}= 'BLUG1 parameter            ';
      phystxt{23}= 'BLUG2 parameter            ';
elseif (iforw==11)
%ORCA
      phystxt{1} = 'Water depth (m)            ';
      phystxt{2} = 'Sediment sound speed (m/s) ';
      phystxt{3} = 'Shear sound speed (m/s)    ';
      phystxt{4} = 'Attenuation (dB/\lambda)    ';
      phystxt{5} = 'Shear Attenuation (dB/\lambda)';
      phystxt{6} = 'Sediment density (g/cm^3)   ';
      phystxt{7} = 'Layer Thickness (m)        ';
      phystxt{8} = 'Source depth (m)           ';
      phystxt{9} = 'Source range (m)           ';
      phystxt{11}= 'Shape coefficient          ';
      phystxt{15}= 'Receiver depth (m)         ';
      phystxt{19}= 'Array tilt (m)             ';
      phystxt{20}= 'Ocean Sound speed  (m/s)   ';
      phystxt{21}= 'Noise                      ';
    end

    xtext=phystxt{iparm};
