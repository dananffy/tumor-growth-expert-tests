
% % 1D random variable following the uniform distribution, data as density
% % convoluted with Gaussians

function [x_out, count] = MetroHast5(y,sigma,nIter)

global N;

x = zeros(nIter,1);
x(1) = 0.7;         
count = 0;

sigma_z = 0.02;     

s1 = [16 20 22 24 24 26 27 28 32]; 
s2 = [20 24 30 26 30 15 20 30 25];

J = length(s1); 
[h1, h2] = meshgrid(-(N-1)/2:(N-1)/2, -(N-1)/2:(N-1)/2);

for k = 2 : nIter
    
      current = x(k - 1);
    
    % Gaussian sampler
      proposed = normrnd(current, sigma_z);  % x'=x+z, where z ~ N(0,sigma_z^2)

      F_p = Test(proposed);  %eta(x') 
      F_c = Test(current);   %eta(x)
    
      F_prop = zeros(N,N);
      F_curr = zeros(N,N);
  
 for i = 1 : J
     hg = exp(-((h1 - s1(i)).^2 + (h2 - s2(i)).^2)/(2*0.1^2));
     F_prop = F_prop + conv2(F_p, hg, 'same');
     F_curr = F_curr + conv2(F_c, hg, 'same');
 end
   
 if 0.5 <= proposed && proposed <= 0.8 && 0.5 <= current && current <= 0.8

     LH1 = exp(-norm(y-F_prop,'fro')^2/2/sigma^2)*0.3;   %pi(x'|y)  
     LH2 = exp(-norm(y-F_curr,'fro')^2/2/sigma^2)*0.3;   %pi(x|y)
      
     A = LH1/LH2;   %y is the observed data vector
     A = min(1,A);
 else 
        A = 0;
 end

    if(rand < A)
         x(k) = proposed;
         count = count + 1;
    else
         x(k) = current;
    end
end

x_out = x(nIter*0.25:nIter);    %discard the burn-in phase






