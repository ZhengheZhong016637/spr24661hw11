function problem1
clear all; close all;
a = sqrt(2);
h = 0.05;
k = h/a;
N = 240;
Nt = 80;
phi = @(x) max(1-abs(x),0);
e = ones(N,1);
o = zeros(N,1);
A = spdiags([e o e],-1:1,N,N);
A(1,N) = 1;
A(N,1) = 1;
B = spdiags([-e o e],-1:1,N,N);
B(1,N) = -1;
B(N,1) = 1;
C = [-sqrt(2),sqrt(2);1,1];
Dp = spdiags([-e e],-1:0,N,N);
Dp(1,N) = -1;
Dn = spdiags([-e e],0:1,N,N);
Dn(N,1) = 1;
E = spdiags([e -2*e e],-1:1,N,N);
E(1,N) = 1;
E(N,1) = 1;
Fp = spdiags([e -4*e 3*e],-2:0,N,N);
Fp(1,N) = -4;
Fp(1,N-1) = 1;
Fp(2,N) = 1;
Fn = spdiags([-3*e 4*e -e],0:2,N,N);
Fn(N,1) = 4;
Fn(N,2) = -1;
Fn(N-1,1) = -1;
Gp = spdiags([e -2*e e],-2:0,N,N);
Gp(1,N-1) = 1;
Gp(1,N) = -2;
Gp(2,N)=1;
Gn = spdiags([e -2*e e],0:2,N,N);
Gn(N-1,1) = 1;
Gn(N,1) = -2;
Gn(N,2) = 1;

solve(4)

    function solve(option)
        xi = initial;
        eta = initial;
        switch option
            case 1
                method = 'Lax Friedrich';
            case 2
                method = 'Upwind';
            case 3
                method = 'Lax Wendroff';
            case 4
                method = 'Beam Warming';            
        end
        for i = 1:Nt
            switch option
                case 1
                    xi = LaxFstep(xi,-a);
                    eta = LaxFstep(eta,a);
                case 2
                    xi = Upwind(xi,-a);
                    eta = Upwind(eta,a);
                case 3
                    xi = LW(xi,-a);
                    eta = LW(eta,a);
                case 4
                    xi = BW(xi,-a);
                    eta = BW(eta,a);
            end
            if i == 10
                u1 = xe2u(xi,eta);
                u1e = exactsol(i*k);
                time = ', t = 1/2a.';
                figure; hold on
                plot(u1e)
                plot(-1*u1)                
                legend('exact','numerical')
                title(strcat(method,time))
            elseif i == 20
                u2 = xe2u(xi,eta);
                u2e = exactsol(i*k);
                time = ', t = 1/a.';
                figure; hold on
                plot(u2e)
                plot(-1*u2)
                legend('exact','numerical')
                title(strcat(method,time))
            elseif i == 40
                u3 = xe2u(xi,eta);
                u3e = exactsol(i*k);
                time = ', t = 2/a.';
                figure; hold on
                plot(u3e)
                plot(-1*u3)
                legend('exact','numerical')
                title(strcat(method,time))
            elseif i == 80
                u4 = xe2u(xi,eta);
                u4e = exactsol(i*k);
                time = ', t = 4/a.';
                figure; hold on
                plot(u4e)
                plot(-1*u4)
                legend('exact','numerical')
                title(strcat(method,time))
            end
        end        
    end
    
    function unew = LaxFstep(u,a)
        unew = 0.5*A*u-(a*k/(2*h))*B*u;
    end

    function unew = Upwind(u,a)             
        if a>0
            unew = u-(a*k/h)*Dp*u;
        else
            unew = u-(a*k/h)*Dn*u;
        end
    end

    function unew = LW(u,a)
        unew = u - (a*k/(2*h))*B*u+(a^2*k^2/(2*h^2))*E*u;
    end

    function unew = BW(u,a)
        if a>0
            unew = u-(a*k/(2*h))*Fp*u+(a^2*k^2/(2*h^2))*Gp*u;
        else
            unew = u-(a*k/(2*h))*Fn*u+(a^2*k^2/(2*h^2))*Gn*u;
        end
    end

    function u = initial
        u = zeros(N,1);
        for i=102:120
            u(i) = -0.5;
        end
        for i = 122:140
            u(i) = 0.5;
        end
    end

    function u = xe2u(xi,eta)
        w = zeros(2,N);
        for j = 1:N
            w(:,j) = C*[xi(j);eta(j)];
        end
        ux = w(2,:);
        u = zeros(N,1);
        for i=2:N
            u(i) = u(i-1)+h*ux(i-1);
        end
    end

    function ue = exactsol(t)
        ue = zeros(N,1);
        for i = 1:N
            ue(i) = (1/2)*(phi((i-1)*h-6+a*t)+phi((i-1)*h-6-a*t));
        end
    end
end