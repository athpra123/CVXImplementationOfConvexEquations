function [B,P] = fn(U,Bmin,Bmax,Pmin,Pmax,t_0,tmin,tmax,alpha,j,N_0,del,d_0)
    Bj = (Bmax - Bmin) * rand(U,1) + Bmin;
    Pj = (Pmax - Pmin) * rand(U,1) + Pmin;
    t = (tmax - tmin) * rand(U, 1) + tmin;
    o = rand(U,1);
    alpha = rand(U,1);
    sgm = 0.5

    cvx_begin
    cvx_solver sedumi
    variable f(U);
    variable y(U);
    variable x(U);
    variable z(U);
    variable m(U);
    variable q(U);
    variable B(U);
    variable P(U);
    
    mj = 5 * rand(U,1);
    xj =  5 * rand(U,1);
    yj = 5 * rand(U,1);
    
    minimize sum(-alpha.*f)
    for i =1:U
        B(i) >= Bmin %% 12a)
        sum(B(i)) <= Bmax %% 12b)
        P(i) >= Pmin %% 12c)
        sum(P(i)) <= Pmax %% 12d)
        
        f(i) <= exp(yj(i)) + (y(i)-yj(i))*exp(yj(i)); %% 26a)
        y(i) <= (-1/2)*(x(i)).^2; %% 26b)
        m(i) >= 2.^q(i)-1 %% 26c)
        q(i) >= d_0*(1-sgm) %% 26d) but removing the division by (B(i)*t(i))
        z(i) >= 1/4*((B(i)+m(i)).^2 - 2*(B(i)-m(i))*(Bj(i)-mj(i))+(Bj(i)-mj(i)).^2) %% 23)
        4*N_0*z(i) <= 2*(x(i)+P(i))*(xj(i)+Pj(i)) - (xj(i)+Pj(i)).^2 - 2*(x(i)-P(i))*(xj(i)-Pj(i)) + (xj(i)+Pj(i)).^2 %% 25)
    end
    cvx_end
end