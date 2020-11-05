function [] = plotHHModel(Vm_vector, n_vector, m_vector, h_vector)
    subplot (2,1,1);
    plot(Vm_vector);

    subplot (2,1,2)
    plot(n_vector);
    hold on;
    plot(m_vector);
    hold on;
    plot(h_vector);
end