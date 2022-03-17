function m = mass_fun(t)
    global m_dot M0 m_fuel;
    
    if m_fuel - (m_dot*t) > 0 
        m =  (m_fuel - (m_dot*t)) + M0;
    else
        m = M0;
    end
    
end