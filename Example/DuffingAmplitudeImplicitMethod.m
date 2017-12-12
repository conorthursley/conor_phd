A=linspace(3,0,300);
w0=1;
h=+1/2;
Fc=1/2;
w=zeros(300,1)';
figure
hold on
for ii=1:300
    a=A(ii);
    w(ii)=abs(((0.353553i)*sqrt(-3*a^3-8*a+4))/sqrt(a));
end

plot(w,A)

