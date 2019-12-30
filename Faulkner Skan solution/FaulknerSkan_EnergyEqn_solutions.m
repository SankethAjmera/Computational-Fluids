clc
clear

% Falker-Skan and energy equation solutions fro wedge flow
% Using 4th order Runge Kutta method and Shooting Technique
% Variable m, Prandtl number, Eckert number and lambda
% Set lambda = 0 for isothermal heating
% Created by Sanketh Ajmera 

% Input example 

% External flow over flat plate (m = 0) for a fluid with Pr = 1 without
% viscous dissipation (Ec = 0) for constant wall temperature (lambda = 0)
m =   0;      % Shape factor, m = 0 for flat plate
Ec =  0;      % Eckert number
Pr =  1;      % Prandtl number
lambda = 0;   % Lambda

% Output example
% eta = eta ; f = f ; G = f' ; H = f" ; I = f"'
% T = Theta ; P = Theta' ; Q = Theta" 

[eta,f,G,H,I] = FaulknerSkan(m); % Soln for Faulker Skan for a given m
[T,P,Q] = EnergyEqn(m,Pr,Ec,lambda); % Soln energy eqn for given m,Pr,Ec,lambda


%%

%----------------------------|Plotting|-----------------------------------%

% Variable m

 m = [0,1/3,2/3,1];
 Pr = 1;
 Ec = 0;
 lambda = 0;

for i = 1:length(m)
    [eta,f,G,H,I] = FaulknerSkan(m(i)); 
    [T,P,Q] = EnergyEqn(m(i),Pr,Ec,lambda);
    plot(eta,G,'DisplayName',['{\it{m}} = ',num2str(round(m(i),2))])
    hold on
    legend;
    xlabel('\eta')
    ylabel('\it{f''}')
    set(gca,'FontSize',12)
    title('f '' vs \eta for external flow over a flat plate without viscous disspiation')
end
hold off

%% Variable Pr

 Pr = [0.1,0.5,1,2,5];
 m = 0;
 Ec = 0;
 lambda = 0;

for i = 1:length(Pr)
    [T,P,Q] = EnergyEqn(m,Pr(i),Ec,lambda);
    plot(eta,Q,'DisplayName',['{\it{Pr}} = ',num2str(round(Pr(i),2))])
    hold on
    legend;
    xlabel('\eta')
    ylabel('\theta2')
    set(gca,'FontSize',12)
    title('Effect of Pr on \theta2')
end
hold off

%% Variable Ec

 Pr = 1;
 m  = 0;
 Ec = [0.1,0.5,1,2,5];
 lambda = 0;

for i = 1:length(Ec)
    [T,P,Q] = EnergyEqn(m,Pr,Ec(i),lambda);
    plot(eta,Q,'DisplayName',['{\it{Ec}} = ',num2str(round(Ec(i),2))])
    hold on
    legend;
    xlabel('\eta')
    ylabel('\theta2')
    set(gca,'FontSize',12)
    title('Effect of Ec on \theta2')
end

%% Variable lambda

 Pr = 1;
 m  = 0;
 Ec = 0;
 lambda = [0,0.5,2,10];

for i = 1:length(lambda)
    [T,P,Q] = EnergyEqn(m,Pr,Ec,lambda(i));
    plot(eta,Q,'DisplayName',['{\lambda} = ',num2str(round(lambda(i),2))])
    hold on
    legend;
    xlabel('\eta')
    ylabel('\theta2')
    set(gca,'FontSize',12)
    title('Effect of \lambda on \theta2')
end

%% All functions defined here

% f' = G                                --(1)
% G' = H                                --(2)
% H' = m(G^2-1)-(m+1)fH/2               --(3)
% T' = P                                --(4)
% P' = Pr*(-Ec*(f")^2 - f/2*(1+m)*P)    --(5)

% Define F1 = G = f'                     =(1)
function [val] = F1(G)
val = G;
end
% Define F2 = H = f''                    =(2)
function [val] = F2(H)
val = H;
end
% Define F3 = H'                         =(3)
function [val] = F3(f,G,H,m)
val = m*(G.^2-1)-(m+1)*f.*H/2;
end
% Define F4 = P = T'                     =(4)
function [val] = F4(P)
val = P;
end
% Define F5 = P' = T"                    =(5)
function [val] = F5(T,P,f,G,H,m,Pr,Ec,lambda)
val = Pr*(-Ec*(H.^2) - f./2*(1+m).*P+lambda*G.*T);
end

% Define function for returning f,f' and f" for a given IC (H(0))
function [x,f,G,H] = freturn(H_Guess,n,h,m)

x = zeros(n,1); x(1) = 0; % Initial condition for eta
f = zeros(n,1); f(1) = 0; % Initial condition for f
G = zeros(n,1); G(1) = 0; % Initial condition for f"
H = zeros(n,1); H(1) = H_Guess; % Guess IC for f"

for i = 2:n
    k1 = h*F1(G(i-1));
    l1 = h*F2(H(i-1));
    m1 = h*F3(f(i-1),G(i-1),H(i-1),m);
    
    k2 = h*F1(G(i-1)+l1/2);
    l2 = h*F2(H(i-1)+m1/2);
    m2 = h*F3(f(i-1)+k1/2,G(i-1)+l1/2,H(i-1)+m1/2,m);
    
    k3 = h*F1(G(i-1)+l2/2);
    l3 = h*F2(H(i-1)+m2/2);
    m3 = h*F3(f(i-1)+k2/2,G(i-1)+l2/2,H(i-1)+m2/2,m);
    
    k4 = h*F1(G(i-1)+l3);
    l4 = h*F2(H(i-1)+m3);
    m4 = h*F3(f(i-1)+k3,G(i-1)+l3,H(i-1)+m3,m);
    
    
    x(i) = x(i-1)+h;                   % eta
    f(i) = f(i-1)+(k1+2*k2+2*k3+k4)/6; % f
    G(i) = G(i-1)+(l1+2*l2+2*l3+l4)/6; % f'
    H(i) = H(i-1)+(m1+2*m2+2*m3+m4)/6; % f"
    
end
end

% Complete solution for faulkner skan equation for any m
function [eta,f,G,H,I] = FaulknerSkan(m)

L = 10;     % eta = 10, same as eta ---> inf
n = 50;     % numer of iterations
h = L/n;    % step size

G0 = 1;     % actual value of G(inf), convergence condition
err = 1e6;  % error
x = 0.1;    % step size for guessing H(0)

H_Guess = 0.5; % Inital guess for H(0)
H1 = 0; % Initial guess for H(0)1
H2 = 2; % Initial guess for H(0)2

[~,~,G_vec1,~] = freturn(H1,n,h,m);
[~,~,G_vec2,~] = freturn(H2,n,h,m);

G1 = G_vec1(end); % G(inf) at H1 = H(0)
G2 = G_vec2(end); % G(inf) at H2 = H(0)

% Attain two reasonable initial guesses for H(0)

while isnan(G1) || G1 < 0
    H1 = H1+x; % increment by x (0.1)
    [~,~,G_vec1,~] = freturn(H1,n,h,m);
    G1 = G_vec1(end);
end
while isnan(G2) || G2 < 0
    H2 = H2-x; % decrement by x (0.1)
    [~,~,G_vec2,~] = freturn(H2,n,h,m);
    G2 = G_vec2(end);
end
H1 = H1+x;
H2 = H2-x;


while err > 1e-6 % Error
    
    [~,~,G_vec1,~] = freturn(H1,n,h,m);
    [~,~,G_vec2,~] = freturn(H2,n,h,m);
    
    G1 = G_vec1(end); % G(inf) at H1 = H(0)
    G2 = G_vec2(end); % G(inf) at H2 = H(0)
    
    H_Guess = (H2-H1)/(G2-G1)*(G0-G1)+H1; % new guess for H(0)
    
    % Bisection method to check where the root lies
    
    if (G2-G0)*(G1-G0)>0
        H1 = H_Guess;
        err = abs(G0-G1);
    else
        H2 = H_Guess;
        err = abs(G0-G2);
    end
    
end


[eta,f,G,H] = freturn(H_Guess,n,h,m);
I = F3(f,G,H,m);

end

% Define function for returning Theta,Theta' for a given IC (Th'(0))
function [T,P] = Treturn(P_guess,f,G,H,n,h,m,Pr,Ec,lambda)

T = zeros(n,1); T(1) = 1;       % Initial condition for theta
P = zeros(n,1); P(1) = P_guess; % Guess IC for theta'

for i = 2:n
    
    k1 = h*F4(P(i-1));
    l1 = h*F5(T(i-1),P(i-1),f(i-1),G(i-1),H(i-1),m,Pr,Ec,lambda);
    
    k2 = h*F4(P(i-1)+l1/2);
    l2 = h*F5(T(i-1)+k1/2,P(i-1)+l1/2,f(i-1),G(i-1),H(i-1),m,Pr,Ec,lambda);
    
    k3 = h*F4(P(i-1)+l2/2);
    l3 = h*F5(T(i-1)+k2/2,P(i-1)+l2/2,f(i-1),G(i-1),H(i-1),m,Pr,Ec,lambda);
    
    k4 = h*F4(P(i-1)+l3);
    l4 = h*F5(T(i-1)+k3,P(i-1)+l3,f(i-1),G(i-1),H(i-1),m,Pr,Ec,lambda);
    
    T(i) = T(i-1)+(k1+2*k2+2*k3+k4)/6; % Theta
    P(i) = P(i-1)+(l1+2*l2+2*l3+l4)/6; % Theta'
    
end
end

% Complete solutin of energy equation for a given m, Pr ,Ec
function [T,P,Q] = EnergyEqn(m,Pr,Ec,lambda)
L = 10;     % eta = 10, same as eta ---> inf
n = 50;     % Numer of iterations
h = L/n;    % step size

[~,f,G,H] = FaulknerSkan(m); % Solve for Faulker Skan for a given m


P_Guess = 1 ;   % initial guess value for Th'(0)
Tinf = 0;       % actual value of T(inf), convergence condition
err2 = 1e6;     % error
y = 0.01;       % increment for guessing appropriate Th'(0) = P(0)
P1 = -2;        % Initial guess lower limit for P(0)1
P2 =  2;        % Initial guess higher limit for P(0)2
T1 = NaN;
T2 = NaN;


while isnan(T1) || abs(T1) < 1e-6
    P1 = P1+y; % increment by x (0.1)
    [T_1,~] = Treturn(P_Guess,f,G,H,n,h,m,Pr,Ec,lambda);
    T1 = T_1(end);
end
while isnan(T2) || abs(T2) < 1e-6
    P2 = P2-y; % increment by x (0.1)
    [T_2,~] = Treturn(P_Guess,f,G,H,n,h,m,Pr,Ec,lambda);
    T2 = T_2(end);
end

P1 = P1+y;
P2 = P2-y;


while err2 > 1e-6 && ~(P1== P2) % Error
    
    [T_vec1,~] = Treturn(P1,f,G,H,n,h,m,Pr,Ec,lambda);
    [T_vec2,~] = Treturn(P2,f,G,H,n,h,m,Pr,Ec,lambda);
    
    T1 = T_vec1(end); % T(inf) at P1 = Th'(0)
    T2 = T_vec2(end); % G(inf) at P2 = Th'(0)
    
    P_Guess = (P2-P1)/(T2-T1)*(Tinf-T1)+P1; % new guess for H(0)
    
    % Bisection method to check where the root lies
    
    if (T2-Tinf)*(T1-Tinf)>0
        P1 = P_Guess;
        err2 = abs(Tinf-T1);
    else
        P2 = P_Guess;
        err2 = abs(Tinf-T2);
    end

end

[T,P] = Treturn(P_Guess,f,G,H,n,h,m,Pr,Ec,lambda); % outputs [Theta,Theta']
Q = F5(T,P,f,G,H,m,Pr,Ec,lambda); % outputs Theta"
end


