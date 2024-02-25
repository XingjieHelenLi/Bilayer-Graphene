function [f e] = Test_Tersoff(x)
    %% usage: 
    % [f e] = Tersoff(x)
    % input, 
    %       x: coordinates atoms, 3 by n
    % output,
    %       f: forces on atoms, 3 by n
    %       e: total potential energy, scalar, a real number
    % 
    % Notice: the forces and total potential energy are calculated analitically
    %         without any boundary conditions.
    %         this matlab script is not supposed to be used for large scale
    %         (more than 10,000 atoms) simulations.
    f = zeros(size(x));
    e = 0;
    
       
    %% parameters for lattice structure 
    % lattice spacing for graphene, change it for other lattice structures 
    a0 = 0.14606009615836723636E1; 
    
    % number of particles
    np = length(x);
    
    % potential energy for each particle, the sum is the total potential
    % energy
    p = ones(np, 1);
    
    %% parameters for tersoff potential (for carbon only)
    % source: J. Tersoff, Phys Rev B, 39, 5566 (1989)
    lam1 = 3.4879;
    lam2 = 2.2119;
    ss   = 2.1;
    a    = 1393.6;
    b    = 346.7;
    mu   = 2.2119;
    beta = 1.5724E-7;
    ci2  = 14.47726401E8;
    di2  = 18.90858256;
    h    = -0.57058;
    betan= 1.12313437E-5;
    pn   = -0.687275827;
    nn   = 0.72751;
    rcut = 2.1;
    rcut2= 4.41;
    
    %% Find the neighboring atoms
    nnbmax = 100;
    mn = zeros(np, 1);
    neib = zeros(np, nnbmax);
    for i = 1:np
        for j = i+1:np
            xij = x(:,i) - x(:,j);
            rij = sqrt(xij(1)*xij(1) + xij(2)*xij(2) + xij(3)*xij(3));
            if(rij < rcut)
                mn(i) = mn(i) + 1;
                if(mn(j) > nnbmax)
                    disp('Please check your input data, too many neighboring atoms');
                    exit;
                end
                neib(i, mn(i)) = j;
                rneib(i, mn(i), 1:3) = xij;
                rneib(i, mn(i), 4) = rij;
                
                mn(j) = mn(j) + 1;
                if(mn(j) > nnbmax)
                    disp('Please check your input data, too many neighboring atoms');
                    exit;
                end
                neib(j, mn(j)) = i;
                rneib(j, mn(j), 1:3) = -xij;
                rneib(j, mn(j), 4) = rij;
            end
        end
    end
    
    %% compute the force and energy
    for i = 1:np
        for ji = 1:mn(i)
            j = neib(i, ji);
            xij(1:3,1) = rneib(i, ji, 1:3);
            rij = rneib(i, ji, 4);
            rij2 = rij*rij;
            urij = xij/rij;
            if(rij < ss)
                fr = a*exp(-lam1*rij); dfr =-lam1*fr;
                fa =-b*exp(-lam2*rij); dfa =-lam2*fa;
                
                fc = fcut(rij, 0);
                dfc= fcut(rij, 1);
                
                % the first term in the potential
                p(i) = p(i) + fc*fr;
                fcr  = dfc*fr + fc*dfr;
                f(:, i) = f(:, i) + fcr*urij;
                f(:, j) = f(:, j) - fcr*urij;
                
                % compute the angle-dependent term
                zeta = 0;
                for ki = 1:mn(i)
                    k = neib(i, ki);
                    if(j ~= k)
                        xik(1:3, 1) = rneib(i, ki, 1:3);
                        rik = rneib(i, ki, 4);
                        urik = xik/rik;
                        if(rik < ss)
                            cijk = dot(urij, urik);
                            trm = di2 + (h-cijk)^2;
                            g = 1 + (ci2/di2) - (ci2/trm);
                            fck = fcut(rik,0);
                            zeta = zeta + fck*g;
                        end
                    end
                end
                
                trm = betan*zeta^nn;
                zetan = 1 + trm;
                bij = zetan^pn;
                
                if(abs(zeta) < 1E-8)
                    dbij = 0;
                else
                    dbij = pn*bij/zetan*nn*trm/zeta;
                end
                
                p(i) = p(i) + bij*fc*fa/2;
                p(j) = p(j) + bij*fc*fa/2;
                
                fca = bij*(fc*dfa + dfc*fa);
                f(:, i) = f(:, i) + fca*urij;
                f(:, j) = f(:, j) - fca*urij;
                
                fcab = fa*fc*dbij;
                
                for ki = 1:mn(i)
                    k = neib(i, ki);
                    if(j ~= k)
                        xik(1:3, 1) = rneib(i, ki, 1:3);
                        rik = rneib(i, ki, 4);
                        rik2= rik*rik;
                        urik = xik/rik;
                        
                        if(rik < ss)
                            xjk = x(:,j) - x(:,k);
                            rjk = sqrt(xjk(1)*xjk(1) + xjk(2)*xjk(2) + xjk(3)*xjk(3));
                            urjk = xjk/rjk;
                            cijk = dot(urij, urik);
                            trm = di2 + (h-cijk)^2;
                            
                            g = 1+(ci2/di2)-(ci2/trm);
                            dg = ci2/trm^2;
                            trm = -2*(h-cijk);
                            dg = dg*trm;
                            
                            fck = fcut(rik, 0);
                            dfck= fcut(rik, 1);
                            
                            s11 = -(fck*dg*cijk)/rij2;
                            s12 = (fck*dg)/(rij*rik);
                            s21 = s12;
                            s22 = (dfck*g*rik - fck*dg*cijk)/rik2;
                            
                            fij = fcab*(s11+s21)*xij;
                            f(:,i) = f(:,i) + fij;
                            f(:,j) = f(:,j) - fij;
                            
                            fik = fcab*(s12+s22)*xik;
                            f(:,i) = f(:,i) + fik;
                            f(:,k) = f(:,k) - fik;
                            
                            fjk =-fcab*s12*xjk;
                            f(:,j) = f(:,j) + fjk;
                            f(:,k) = f(:,k) - fjk;
                        end
                    end
                end
                            
            end
        end
    end
    
    e = sum(p);
    
end

%% local cut-off function in tersoff potential
function fc = fcut(r, iflag)
    
    % parameters for tersoff potential 
    rr   = 1.8;
    ss   = 2.1;
    
    if(iflag == 0)
        if(r < rr)
            fc = 1;
        elseif(r < ss)
            trm = pi*(r-rr)/(ss-rr);
            fc = 0.5 + 0.5*cos(trm);
        else
            fc = 0;
        end
    elseif(iflag == 1)
        if(r < rr)
            fc = 0;
        elseif(r < ss)
            trm = pi*(r-rr)/(ss-rr);
            fc =-0.5*pi/(ss-rr)*sin(trm);
        else
            fc = 0;
        end
    end
    
end
