%% Test 1(b) in the paper 

clear all;
clc;

sigma = 0.125;
nIter = 800;
global mu c0 Final_Time N;

%% set up parameters   

mu = 0.5;         %mean of prior N(mu,c0^2)
c0 = 0.5;         %std of prior

G_star = 1;       %underlying true G
Final_Time = 0.5; 

dx = 0.1;
a = -2.2; b = 2.2;                 %default domain
[X,Y] = meshgrid(a:dx:b,a:dx:b);
xx = [a:dx:b];

%% create observed data
s1 = [16 20 22 24 24 26 27 28 32]; 
s2 = [20 24 30 26 30 15 20 30 25];

J = length(s1); 
y0 = Test(G_star);
N = size(y0,1);
y = zeros(size(y0));
[h1, h2] = meshgrid(-(N-1)/2:(N-1)/2, -(N-1)/2:(N-1)/2);

for i = 1 : J
   hg = exp(-((h1 - s1(i)).^2 + (h2 - s2(i)).^2)/(2*0.1^2));    
   y = y + conv2(y0, hg, 'same'); 
end

P = zeros(nIter,1);   
[P, count] = MetroHast5(y,sigma,nIter);
MSE = sum((P - G_star).^2)/length(P)




