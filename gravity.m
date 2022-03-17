function FG = gravity(r,m)
    global mu re J;
    
    g = [r(1)*(1+((3/2)*J*(re/norm(r))^2))*(1-5*(r(3)/norm(r))^2) ;
        r(2)*(1+((3/2)*J*(re/norm(r))^2))*(1-5*(r(3)/norm(r))^2);
        r(3)*(1+((3/2)*J*(re/norm(r))^2))*(3-5*(r(3)/norm(r))^2)] ;
    FG = -mu*m*(1/(norm(r)^3)).*g;
end

