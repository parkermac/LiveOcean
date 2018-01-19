function [nclongname,ncunits,nctimename] = atm_attributes(var);
% atm_attributes.m  5/9/2014  Parker MacCready
%
% used by make_atm_main.m

switch var
    case 'Pair'
        nclongname = 'surface air pressure';
        ncunits = 'millibar';
        nctimename = 'pair_time';
        scalefactor = 1/100; % convert Pa to mbar
    case 'rain'
        nclongname = 'rain fall rate';
        ncunits = 'kilograms meter-2 second-2';
        nctimename = 'rain_time';
    case 'swrad'
        nclongname = 'solar shortwave radiation flux';
        ncunits = 'watts meter-2';
        nctimename = 'srf_time';
    case 'lwrad_down'
        nclongname = 'downwelling longwave radiation flux';
        ncunits = 'watts meter-2';
        nctimename = 'lrf_time';
    case 'Tair'
        nclongname = 'surface air temperature';
        ncunits = 'Celsius';
        nctimename = 'tair_time';
    case 'Qair'
        nclongname = 'surface air relative humidity';
        ncunits = 'percentage';
        nctimename = 'qair_time';
    case 'Uwind'
        nclongname = 'surface u-wind component';
        ncunits = 'meter second-1';
        nctimename = 'wind_time';
    case 'Vwind'
        nclongname = 'surface v-wind component';
        ncunits = 'meter second-1';
        nctimename = 'wind_time';
end