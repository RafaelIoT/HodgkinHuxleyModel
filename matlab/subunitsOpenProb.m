% Probability of opening channels subunits for the given transition rates
% Inputs: 
%    transition rates
% Outputs: probs - Vector with probabilities of opening channels subunits
% Date: 9 out 2020
% Authors:
%   Rafael Cruz, 50380
%   Diana Castaneda, 51549
function [probs] = subunitsOpenProb(an, bn, am, bm, ah, bh)
    n = an / (an + bn);
    m = am / (am + bm);
    h = ah / (ah + bh);
    
    probs = [n, m, h];
end