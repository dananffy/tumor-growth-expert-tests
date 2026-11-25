%The test when G(x) is a random field

function rho = Random_Field(gf)

global Final_Time a b C1 C2;

initialFun=@(X,Y)(0.9*(((X-C1).^2+(Y-C2).^2)<0.3));

[X,Y,rho,p] = solver_2D('T', Final_Time, 'GrowthFun', gf, 'InitialFun', initialFun, 'a', a, 'b', b);

    