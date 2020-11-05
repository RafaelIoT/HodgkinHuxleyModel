% hodgkinHuxleyModel(6.3, 53, 2, 0.2, 15, 0.01)
function [vectors] = hodgkinHuxleyModel(T, Is, Is_begin, Is_duration, total_time, step)
    Vr = -60;
    Vm0 = -60;
    
    Na_in_concentration = 49.5; % mM;
    Na_ex_concentration = 437.0; % mM;
    K_in_concentration = 397; % mM;
    K_ex_concentration = 29; % mM;

    
    EK = nernstPotential(K_in_concentration, K_ex_concentration, 1, T, "C");
    ENa = nernstPotential(Na_in_concentration, Na_ex_concentration, 1, T, "C");
    
    % EK = -72.100; % mV
    % ENa = 52.4;
    EL = -49.187;

    Cm = 1.0; % micro Farad / cm^2;

    % T = 6.3;
    
    n0 = 0.31768;
    m0 = 0.05293;
    h0 = 0.59612;

    gK_max = 36.0; % mS / cm^2;
    gNa_max = 120.0;
    gL = 0.3;
    gK = 0.367;
    gNa = 0.011;

    % R = 8.3145; % constante dos gases perfeitos J * K^⁻1 * mol^-1;
    % F = 96485; % 96485.332 constante de faraday C / mol;

    time = total_time; % 10 msec 
    delta_t = step; % msec
    total_time_deltas = time / delta_t;
    
    % Is_begin = 2;
    % Is = 53;
    % Is_duration = 0.2;

    Im_vector = stimulusPlot(Is, Is_begin, Is_duration, time, delta_t);

    n_vector = zeros(1, total_time_deltas+1);
    m_vector = zeros(1, total_time_deltas+1);
    h_vector = zeros(1, total_time_deltas+1);

    IK_vector = zeros(1, total_time_deltas+1);
    INa_vector = zeros(1, total_time_deltas+1);
    IL_vector = zeros(1, total_time_deltas+1);

    Vm_vector = zeros(1, total_time_deltas+1);

    gK_vector = zeros(1, total_time_deltas+1);
    gNa_vector = zeros(1, total_time_deltas+1);


    initialRates = channelRates(0);
    initialRatesCell = num2cell(initialRates);
    [an, bn, am, bm, ah, bh] = initialRatesCell{:};

    initialProbs = subunitsOpenProb(an, bn, am, bm, ah, bh);
    initialProbsCell = num2cell(initialProbs);
    [n, m, h] = initialProbsCell{:};

    n_vector (1) = n;
    m_vector (1) = m;
    h_vector (1) = h;

    INa = gNa * (Vr - ENa);
    IK = gK * (Vr - EK);
    IL = gL * (Vr - EL);

    IK_vector (1) = IK;
    INa_vector (1) = INa;
    IL_vector (1) = IL;

    Vm_vector (1) = Vm0;

    gK_vector (1) = gK;
    gNa_vector (1) = gNa;


    for t=1 : total_time_deltas
        vM = Vm_vector(t) - Vr;

        rates = channelRates(vM);
        ratesCell = num2cell(rates);
        [an, bn, am, bm, ah, bh] = ratesCell{:};

        gK = gK_max * n_vector(t)^4;
        gNa = gNa_max * h_vector(t) * m_vector(t)^3;

        INa = gNa * (Vm_vector(t) - ENa);
        IK = gK * (Vm_vector(t) - EK);
        IL = gL * (Vm_vector(t) - EL);

        I = Im_vector(t) - IK - INa - IL;

        Vm_vector(t+1) = Vm_vector(t) + delta_t * I / Cm;
        n_vector(t+1) = n_vector(t) + delta_t * (an * (1 - n_vector(t)) - bn * n_vector(t));
        m_vector(t+1) = m_vector(t) + delta_t * (am * (1 - m_vector(t)) - bm * m_vector(t));
        h_vector(t+1) = h_vector(t) + delta_t * (ah * (1 - h_vector(t)) - bh * h_vector(t));

        IK_vector(t) = IK;
        INa_vector(t) = INa;
        IL_vector(t) = IL;
        
        gK_vector(t) = gK;
        gNa_vector(t) = gNa;
    end
    
    vectors.gK = gK_vector;
    vectors.gNa = gNa_vector;
    vectors.IK = IK_vector;
    vectors.INa = INa_vector;
    vectors.IL = IL_vector;
    vectors.Vm = Vm_vector;
    vectors.Im = Im_vector;
    vectors.n = n_vector;
    vectors.m = m_vector;
    vectors.h = h_vector;
    
end
