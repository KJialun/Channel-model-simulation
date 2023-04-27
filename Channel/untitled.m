clear 

x = [-3:.1:3];
y = zeros(1,length(x));
z = normpdf(x,0,1);
figure
plot3(x,y,z)
hold on
y = [-3:.1:3];
x = zeros(1,length(y));
plot3(x,y,z)
xlabel('X ')
ylabel('Y')
zlabel('PDF')

r=sqrt(x.^2+y.^2);
plot3(x,y,z)