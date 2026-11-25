%% Test 3 in the paper when G(x) is a random field

close all;
clc;

global d sigma nIter Final_Time;
global mu c1 c2 c3;
global X Y a b C1 C2 G0;

d = 3;

%choice = 1;     %rho no noise
choice = 2;      %rho perturbed with Gaussian noise
%choice = 3;      %rho convoluted with Gaussians

%% set up parameters  

mu = 0;          %mean of prior N(mu,c_i^2)
c1 = 0.4;       %std of prior
c2 = 0.3;
c3 = 0.2;
G0 = 2;

sigma = 0.125;

C1 = 1;  %center of the initial circle 
C2 = 1;  

nIter = 500;
Final_Time = 0.5;       

a = -0.5; b = 2.5;
dx = 0.1;           
[X,Y] = meshgrid(a:dx:b,a:dx:b);
xx = [a:dx:b];

M = 2;
s = 2;
cc = 4*pi^2/M^2;

Gamma_1 = (cc*2)^(-s/2);  %i_1=1, i_2=1
Gamma_2 = (cc*5)^(-s/2);  %i_1=1, i_2=2
Gamma_3 = (cc*5)^(-s/2);  %i_1=2, i_2=1
Gamma = [Gamma_1 Gamma_2 Gamma_3];

Zeta_star = [0.5 0.3 0.2];
g_star = Zeta_star.*Gamma;     %true data for g
G_star = G_Function(g_star);   %true G function

%% create observed data

y0 = Random_Field(G_star); 
N = size(y0,1);

if choice == 1
    y = y0;                       %Case I: generate data with true G_star, without noise
elseif choice == 2
    y = y0 + sigma*randn(N,N);    %Case II: generate data with true data and add a Gaussian noise N(0,sigma^2)
else
    y = zeros(N,N);
 for i = 1 : J
    hg = exp(-((h1-s1(i)).^2 + (h2-s2(i)).^2)/(2*0.25^2));
    y = y + conv2(y0,hg,'same');
 end
    y = norm(y,1);  
end

[P, count] = MetroHast_RandomField(y); 

pp = length(P(:,1));
g_1 = mean(P(:,1))
g_2 = mean(P(:,2))
g_3 = mean(P(:,3))
g = [g_1, g_2, g_3];

G1 = G_star(0,X,Y);    %true G
G2 = G_Fun(g);         %Posterior mean of G

GG = zeros(N,N,pp); 
G_temp = zeros(N,N);
E = zeros(pp,1);

for j = 1 : pp      
   E(j) = norm(G_Fun(P(j,:))- G1,'fro')/norm(G1,'fro'); 
end

MSE = sum(E.^2)/pp












