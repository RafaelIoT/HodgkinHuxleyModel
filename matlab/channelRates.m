function [rates] = channelRates(vM)

    if vM == 10 
        an = 0.1;
    else
        an = .01 * ((10-vM) / (exp((10-vM)/10)-1));
    end
    
    bn = .125*exp(-vM/80);
    
    if vM == 25
        am = 0.0;
    else
        am = .1*((25-vM) / (exp((25-vM)/10)-1));
    end
    bm = 4*exp(-vM/18);

    ah = 0.07*exp(-vM/20);
    bh = 1/(exp((30-vM)/10)+1);
    
    rates = [an, bn, am, bm, ah, bh];
end