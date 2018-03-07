%% AMM spring matrix allocator
function [K] = spring_matrix(n,k1,k2)

K=zeros(n*2);
K(1,1)=k2+(k1);
K(1,2)=-k2;
K(2,1)=-k2;
K(2,2)=k2;
%---------------------------------------------------
% last 2x2 diagonal, uses the last cell's secondary mass, y(end)
% populate
if n>1
    q=2*n;
    K(q,q)=k2;
    K(q-1,q-1)=k2+k1;
    K(q-1,q)=-k2;
    K(q,q-1)=-k2;
    K(q-1,q-2)=-k1;
    K(q-2,q-1)=-k1;
end
%---------------------------------------------------
% Iteration of remainder of matrix, ignoring the first row and last diagonal
if n>2
    
    for i=2:((2*n)-2)
        if mod(i,2)==1  % odd number and therefore a primary mass
            K(i,i)=k2+2*k1;
            K(i,i-1)=-k1;
            K(i-1,i)=-k1;
        elseif mod(i,2)==0 % even number and therefore a secondary mass
            K(i,i)=k2;
            K(i,i-1)=-k2;
            K(i-1,i)=-k2;
        end
    end
end
%---------------------------------------------------
end
