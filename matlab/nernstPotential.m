% Nernst potential 
% Inputs: 
%    T - temperature
%    temp_unit - "C"/"K"/"F"
%    q_in_conc - Intracelullar Ion Concentration, mM
%    q_ex_conc - Extracelullar Ion Concentration, mM
% Outputs: E - Nernst potential in mV
% Date: 9 out 2020
% Authors:
%   Rafael Cruz, 50380
%   Diana Castaneda, ....

function [E] = nernstPotential(q_in_conc, q_ex_conc, z, temp, temp_unit)
    % Example: nernstPotential(0.22, 0.02, 1, 26.85, "C")
    
    R = 8.3145; % constante dos gases perfeitos J * K^⁻1 * mol^-1;
    F = 96485; % 96485.332 constante de faraday C / mol;
    if strcmp(temp_unit,"K") == 0
        if strcmp(temp_unit,"C")
            temp = temp + 273.15;
        elseif strcmp(temp_unit,"F")
            temp = 273.5 + ((temp - 32.0) * (5.0/9.0));
        else
            disp("Wrong temperature unit: K, C or F")
        end
    end
    div_conc = q_in_conc / q_ex_conc;
    E = -1000*R*temp/(z*F) * log(div_conc);
end