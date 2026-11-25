%The test when G(x) is a Gaussian field 

function [x_out, count] = MetroHast_RandomField(y)

global d sigma nIter;
global mu c1 c2 c3;

x = zeros(nIter,d);  
x(1,1) = 0.02;
x(1,2) = 0.01;
x(1,3) = 0.01;

count = 0; 
sigma_z = 0.005;

for k = 2 : nIter
    
    current = x(k-1,:);  
    proposed = normrnd(current, sigma_z);  %coefficient vector Xi
    
    G_proposed = G_Function(proposed);   %return G(x)
    G_current = G_Function(current);
    
    F_prop = Random_Field(G_proposed);   %return rho
    F_curr = Random_Field(G_current);      

    LH1 = exp(-norm(y-F_prop,'fro')^2/2/sigma^2)*exp(-(proposed(1)-mu)^2/2/c1^2-(proposed(2)-mu)^2/2/c2^2-(proposed(3)-mu)^2/2/c3^2);
    LH2 = exp(-norm(y-F_curr,'fro')^2/2/sigma^2)*exp(-(current(1)-mu)^2/2/c1^2-(current(2)-mu)^2/2/c2^2-(current(3)-mu)^2/2/c3^2);

    A = LH1/LH2;       
    A = min(1,A);
  
   if(rand < A)
        x(k,:) = proposed;
        count = count + 1;
    else
        x(k,:) = current;
    end
end
    
    x_out = x(ceil(nIter*0.25):nIter,:);  




