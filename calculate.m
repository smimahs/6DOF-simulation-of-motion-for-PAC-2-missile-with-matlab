function y_dot = calculate(t,y)
    global w Thrust e1 e2 J1 X_BC J2 l d i acceleration_g acceleration_nong Ixx Iyy Izz m acceleration_body m_dot FA MA ...
        Tempreture a_atmosisa P rho M0 t_end_thrust i_end;
    
    r = [y(1) y(2) y(3)]; %%position in body
    
    q=quatnormalize(y(10:13,:)');
    C_b2I = quat2dcm([q(1) q(2) q(3) q(4)]);
    C_I2B = quat2dcm([-q(1) q(2) q(3) q(4)]);
    
    r_eci = eci2lla(r,[2021 1 1 12 00 00]);
    h = r_eci(3);
    [Tempreture(i,1), a_atmosisa(i,1), P(i,1), rho(i,1)] = atmosisa(h);
    Tempreture(i,2) = h;
    a_atmosisa(i,2) = h;
    P(i,2) = h;
    rho(i,2) = h;
    
    m(i,1) = mass_fun(t); % mass
    m(i,2) = t;
        
    
    V = y(4:6); % velocity in body
    y_dot(1) = V(1);
    y_dot(2) = V(2);
    y_dot(3) = V(3);

    ve = 3*1000; %m/s
    pe = 1.0132533539 * 0.8; % bar
    Ae = (0.87^2)/2;
    FG = gravity(r,m(i,1));
    Thrust(i,1)= -((pe - P(i,1)/100000)*Ae - (m_dot*ve));
    if m(i,1) == M0
        %FG = -gravity(r,m(i,1));
        Thrust(i,1)=0;
        if t_end_thrust == -1 || i_end == -1
            t_end_thrust = t;
            i_end = i;
        end
    end

    Thrust(i,2)= t;

    Cd = 0.26; % alpha = 1 and mach = 3
    CL = 0.12; % alpha = 1 and mach = 3
    Cy = -0.16; % alpha = 1 and mach = 3
    A = 0.87^2;
    %rho(i,1) = 1.225;
    FD =  0.5 .* rho(i,1) .* (norm(V./a_atmosisa(i,1))^2) .* Cd .* A;
    FL = 0.5 .* rho(i,1) .* (norm(V./a_atmosisa(i,1))^2) .* CL .* A;
    FY = 0.5 .* rho(i,1) .* (norm(V./a_atmosisa(i,1))^2) .* Cy .* A;
    
    eulZYX = quat2eul(q);
    
    FA(i,1:3) = [FL;FY;FD];
    FA(i,4) = t;
    MA(i,1:3)= X_BC .*FA(i,1:3);
    MA(i,4) = t;
    
    xe = 2.5146;
    ye = 0;
    ze = 0;
    
    w = [y(7) y(8) y(9)];
    
    Ixx(i,1) = ((1/4)*m(i,1)*(d/2)^2)+((1/12)*m(i,1)*(l)^2);
    Ixx(i,2) = t;
    Iyy(i,1) = ((1/4)*m(i,1)*(d/2)^2)+((1/12)*m(i,1)*(l)^2);
    Iyy(i,2) = t;
    Izz(i,1)= ((1/2)*m(i,1)*(d/2)^2);
    Izz(i,2) = t;
    
    Roughness_angle = [cos(e2)*cos(e1);
                       cos(e2)*sin(e1);
                       sin(e2)];
  
    temp1 = [Ixx(i,1);
            Iyy(i,1) - (m(i,1)*X_BC^2);
            Izz(i,1) - (m(i,1)*X_BC^2)];
    
    temp2 = [w(2)*w(3)*(-Izz(i,1)+Iyy(i,1));
            w(1)*w(3)*(Izz(i,1)-Ixx(i,1));
            w(2)*w(1)*(Ixx(i,1)-Iyy(i,1))];
    
    temp3 = [ye*sin(e2) - ze*cos(e2)*sin(e1);
            ze*cos(e2)*cos(e1) - xe*sin(e2);
            xe*cos(e2)*sin(e1)-ye*cos(e1)*cos(e2)];
        
    temp4 = [0;w(2);w(3)];
    
    temp5 = [0 0 0;
            0 0 -X_BC;
            0 X_BC 0];
    
    temp6 = [w(3)^2 + w(2)^2;
               -(w(1)*w(2));
               -(w(1)*w(3))];
           
    temp7 = [0;w(3);-w(2)];

    w_dot = (MA(i,1:3)' + temp2 + (Thrust(i,1) .* temp3) - (J2 .* temp4) + (temp5 *(FA(i,1:3)' + (Thrust(i,1).*Roughness_angle))) + (m(i,1)*X_BC.*temp6) - (J1.*temp7)) ./ temp1;

    
    matrix1 = [w(3)^2 + w(2)^2;
               - (w_dot(3) + (w(1)*w(2)));
               w_dot(2) - (w(1)*w(3))];
   
    matrix2 = [0;w(3);-w(2)];

    
    
    a = ((1/m(i,1)) .*  FG)-[0;0;1] + (C_b2I *(((1/m(i,1)).*FA(i,1:3)') + ((Thrust(i,1)/m(i,1)).*Roughness_angle) + (X_BC .* matrix1) - ((J1/m(i,1)).*matrix2 )));

    acceleration_g(i,1:3) = (FG ./ m(i,1))-[0;0;1];
    acceleration_g(i,4) = t;
    acceleration_nong(i,1:3) = C_I2B*(C_b2I *(((1/m(i,1)).*FA(i,1:3)') + ((Thrust(i,1)/m(i,1)).*Roughness_angle) + (X_BC .* matrix1) - ((J1/m(i,1)).*matrix2 )));
    acceleration_nong(i,4) = t;
    acceleration_body(i,1:3) = C_I2B * a;
    acceleration_body(i,4) = t;
    
    y_dot(4) = a(1);
    y_dot(5) = a(2);
    y_dot(6) = a(3);
        
    y_dot(7) = w_dot(1);
    y_dot(8) = w_dot(2);
    y_dot(9) = w_dot(3);
       
    y_dot(10)=(-w(1)*y(11) - w(2)*y(12) - w(3)*y(13))/2;
    y_dot(11)=(w(1)*y(10) + w(3)*y(12) - w(2)*y(13))/2;
    y_dot(12)=(w(2)*y(10) - w(3)*y(11) + w(1)*y(13))/2;
    y_dot(13)=(w(3)*y(10) + w(2)*y(11) - w(1)*y(12))/2;
    y_dot = y_dot';

    i = i+1;

end
