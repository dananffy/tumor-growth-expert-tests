%Given the coefficient vector g, compute G(x)

function G_out = G_Function(g)

global X Y G0; 

M = 2;
C = 2*pi/M; 
  
Phi_1 = sin(C*X).*sin(C*Y);     %i_1=1, i_2=1
Phi_2 = sin(C*X).*sin(2*C*Y);   %i_1=1, i_2=2
Phi_3 = sin(2*C*X).*sin(C*Y);   %i_1=2, i_2=1

G_out = @(p,X,Y)(0*p+0*X.^2+G0+...
                     g(1)*Phi_1+g(2)*Phi_2+g(3)*Phi_3); 



