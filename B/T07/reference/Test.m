% % For Test 1

function rho = Test(xi)    

global Final_Time;

initialFun=@(X,Y)(0.9*((sqrt(X.^2+Y.^2)-0.5-sin(4*atan(Y./X))/2)<0));
gf=@(p,X,Y)(0*p+xi);

[X,Y,rho,p] = solver_2D('T', Final_Time, 'GrowthFun', gf, 'InitialFun', initialFun);
