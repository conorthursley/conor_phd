function my_phase()

IC = [1 1; 1 2;1 3;1 4];
hold on
for ii = 1:length(IC(:,1))
    [~,X] = ode45(@EOM,[0 50],IC(ii,:));
    u = X(:,1);
    w = X(:,2);
    plot(u,w,'r')
end
xlabel('u')
ylabel('w')
grid
end

function dX = EOM(t, y)
dX = zeros(2,1);
u  = y(1);
w  = y(2);
A  = 1;
B  = 1;
dX = [w*u^2 - B*u;...
      A - w - w*u^2];
end
