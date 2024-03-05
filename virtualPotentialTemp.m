function thetav = virtualPotentialTemp(T, Td, P)
%{
This function calculates the virtual potential temperature.
INPUTS:
T  =  temperature (K/C)
Td =  dewpoint (K/C)
P  =  pressure (Pa)
%}

% Check for Celsius
if nanmean(T) < 100
    T = T+273.15;
end
if nanmean(Td) < 100
    Td = Td+273.15;
end

P0 = 101325; % reference pressure

% Calculate vapor pressure - see https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
Pvap = 6.11 * 10.^(7.5*(Td - 273.15)./(237.3 + (Td - 273.15))); % outputs Pvap in mbar

% Calculate the mixing ratio in g/kg - see https://www.weather.gov/media/epz/wxcalc/mixingRatio.pdf
r = 621.97 * Pvap ./ (P/100 - Pvap); % /100 to convert P from Pa to mbar

% Calculate potential temperature
theta = T .* (P0 ./ P) .^ 0.286;
thetav = theta .* (1 + 0.61*r/1000); % convert r from g/kg to kg/kg