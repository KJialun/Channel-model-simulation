function [outputArg1,outputArg2] = circle(A,origin)


% 假定要圆心为(30,40),半径为10
[x,y]=meshgrid((1:100)-40,(1:100)-30);
B=sqrt(x.^2+y.^2)<=10;
C=A(B);


end
x = -5:1:5;
y = -5:1:5;
[X Y] = meshgrid(x,y);
r = 10; % radius
indicator = X.^2 + Y.^2 < r^2 ;
X(indicator)=0;