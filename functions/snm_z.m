function m = snm_z(x,dim)
    x_z = atanh(x);
    m_z = snm(x_z,dim);
    m = tanh(m_z);
end

% do median in z-transform space ?