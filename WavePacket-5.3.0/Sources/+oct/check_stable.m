%--------------------------------------------------------------------------
%
% Check stability of (sparse) system matrix A: 
% real part of all eigenvalues should be negative
%
%--------------------------------------------------------------------------

function check_stable ( label, A, N )

% Construct matrix and calculate eigenvalues
matrix = A; % + A';
if nargin>2
    for d=1:length(N)
        matrix = matrix + N{d} * N{d}';
    end
end
eigval=sort(real(eig(full(matrix))),'descend');

% Check for non-negative real part
epsilon = eps('single');
nnrp = nnz(eigval>=-epsilon);
util.disp([num2str(nnrp) ' non-negative eigenvalues found for ' label ])

% Display eigenvalue with largest real part
if nnrp
    if nnrp < 10
        util.disp('Eigenvalue(s) with largest real part(s)')
        for k=1:nnrp
            util.disp ([int2str(k) ' : ' num2str(eigval(k))])
        end
    else
        for k=1:5
            util.disp ([int2str(k) ' : ' num2str(eigval(k))])
        end
        util.disp('... ...')
        for k=nnrp-4:nnrp
            util.disp ([int2str(k) ' : ' num2str(eigval(k))])
        end
        
    end
end

end

