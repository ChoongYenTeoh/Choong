
clc;

N   = 10;       % the matrix is N x N
% r   = 2;        % the rank of the matrix
%df  = 0:0.1:0.9;
df  = linspace(0,0.9,10)
X = randi([0 1],N)

for i = 1:length(df)
%nSamples    = 2*df; % number of observed entries
nSamples = (N^2)-(df(i).*(N^2))

% For this demo, we will use a matrix with integer entries
% because it will make displaying the matrix easier.

rPerm   = randperm(N^2) % use "randsample" if you have the stats toolbox

%random samples removed
% omega = sort(rPerm(1:nSamples));
%specific columns removed
omega = sort(rPerm);
omega = omega(1:nSamples);

Y = NaN(N)
Y(omega) = X(omega)


observations = X(omega)    % the observed entries
mu           = 1;        % smoothing parameter

% The solver runs in seconds
tic
Xk = solver_sNuclearBP( {N,N,omega}, observations, mu)
toc

sprintf('%10.1f',Xk')
fprintf('Relative error, no rounding: %.8f%%\n', norm(X-Xk,'fro')/norm(X,'fro')*100 )
disp( round(Xk*10000)/10000 )
H9 (i)= norm(X-Xk,'fro')/norm(X,'fro')*100
NoC9 (i) = nSamples ./ 10

end
% disp('The "NaN" entries represent unobserved values');
% disp(Y)
% disp('Recovered matrix (rounding to nearest .0001):')
% disp( round(Xk*10000)/10000 )
% % and for reference, here is the original matrix:
% disp('Original matrix:')
% disp( X )

% The relative error (without the rounding) is quite low:
% fprintf('Relative error, no rounding: %.8f%%\n', norm(X-Xk,'fro')/norm(X,'fro')*100 );
% 
% H = norm(X-Xk,'fro')/norm(X,'fro')*100
% NoC = nSamples ./ 10

plot(NoC9,H9)
xlabel('Number of Observable Samples');
ylabel('Estimation Error (%)');
title('Error between estimate and known solution');
xlim([0 10]);
ylim([0 100]);
grid;