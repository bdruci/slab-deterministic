%Energy Independent, 1-D slab with monodirectional beam, multiple regions
%Isotropic scattering (will upload more general case w/ legendre moments)
%REQUIRES: lgwt.m (gets weights and angles)

%INPUT
X = 10; %thickness (cm)
N_x = 1000; %# of x divisions
Sigt = 1*ones(N_x,1); %define Sigt, Sigs for the regions
Sigs = ones(N_x,1);
Sigs(1:N_x/2) = Sigs(1:N_x/2)*0.5; %ex: Sigma_scatter is .5 in the first half, 
Sigs(N_x/2+1:N_x) = Sigs(N_x/2+1:N_x)*0.5;%but .7 in the second region. Same Sig_total
N_mu = 16; %# of ordinates
Q = zeros(N_x,2*N_mu); %external source
e = .00001; %convergence requirement (not implemented)

%create and initialize vectors
dx = X/N_x;


[mu_all, w_all] = lgwt(2*N_mu,-1,1);
mu = mu_all(1:N_mu)';
w  = w_all(1:N_mu)';

psi_pos = zeros(N_x, N_mu);
psi_neg = zeros(N_x, N_mu);
psi_half = zeros(N_x+1,2*N_mu);
q = zeros(N_x,2*N_mu);
q_old = ones(N_x,2*N_mu);

%Bndry Conds
%monodirectional beam
for n = 1:N_mu
    psi_half(1,n) = 0;
    psi_half(N_x+1,N_mu+n) = 0;
end
psi_half(1,1) = 2*N_mu;

c = 0;
psi_pos_old = zeros(N_x,N_mu);

%Srce
%Q = dx*ones(N_x,N_mu)

%begin alg
while(c < 10000) %any other suggestions?
    %forward
    for i = 1:N_x
        for m = 1:N_mu
            top = psi_half(i,m) + q(i,m)*dx./(2*mu(m));
            bot = 1 + Sigt(i)*dx./(2*mu(m));
            psi_pos(i,m) = top./bot;
            psi_half(i+1,m) = 2*psi_pos(i,m)- psi_half(i,m);
        end
    end
    
    %backwards
    for i = N_x:-1:1
        for m = 1:N_mu
            top = psi_half(i+1,N_mu + m) + q(i,N_mu+m)*dx./(2*mu(m));
            bot = 1 + Sigt(i)*dx./(2*mu(m));
            psi_neg(i,m) = top./bot;
            psi_half(i,N_mu + m) = 2*psi_neg(i,m)- psi_half(i+1,N_mu+m);
        end
    end
    
    %scatter source
    qs = zeros(N_x,2*N_mu);
    for i = 1:N_x
        for m = 1:N_mu
            Psi = 0;
            for n = 1:N_mu %gaussian quadrature
                Psi = Psi + w(n)*psi_pos(i,n);
                Psi = Psi + w(n)*psi_neg(i,n);
            end
            qs(i,m) = Psi*Sigs(i)*0.5;
            qs(i,N_mu+m) = Psi*Sigs(i)*0.5;
        end
    end
    
    %update source
    q = qs + Q;
    c = c + 1;
end


%calculate scalar flux
phi = zeros(N_x+1,1);
leak = 0;
inci = 0;
for i = 1:(N_x+1)
    for m = 1:N_mu
        phi(i) = phi(i) + w(m)*(psi_half(i,m) + psi_half(i,m+N_mu));
        if(i == N_x+1)
            leak = leak + mu(m)*w(m)*psi_half(i,m);
        end
        if(i == 1)
            inci = inci + mu(m)*w(m)*psi_half(i,m);
        end
    end
end

xp = 0:dx:X;
plot(xp, phi/(phi(1)));

leak/inci