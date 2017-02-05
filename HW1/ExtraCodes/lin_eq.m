function x = lin_eq(A,B)
%This function finds the solution to Ax = B
[M,N] = size(A);
if size(B,1) ~= M
error('Incompatible dimension of A and B in lin_eq()!')
end
if M == N, x = A^-1*B; %x = inv(A)*B or gaussj(A,B); %Eq.(2.1.1)
elseif M < N %Minimum-norm solution (2.1.7)
x = pinv(A)*B; %A’*(A*A’)^-1*B; or eye(size(A,2))/A*B
else %LSE solution (2.1.10) for M > N
x = pinv(A)*B; %(A’*A)^-1*A’*B or x = A\B
end