function [probs] = subunitsOpenProb(an, bn, am, bm, ah, bh)
    n = an / (an + bn);
    m = am / (am + bm);
    h = ah / (ah + bh);
    
    probs = [n, m, h];
end