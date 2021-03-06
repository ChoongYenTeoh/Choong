close all;
clear;
clc;

% the matrix is N x N
N   = 500;       


df  = linspace(0,0.9,10)

%number of experiments desired to be run
Number_of_Experiments = 10

%this loop restarts the process with a new matrix and saves all the data at the same time.
for j = 1:Number_of_Experiments 
%  matrix generated
    X = randi([0 1],N) 

    
%this loop ensures and saves all the data of the matrix samples/columns being removed    
for i = 1:length(df) 
% number of observed entries
nSamples = (N^2)-(df(i).*(N^2))

% use "randsample" if you have the stats toolbox
rPerm   = randperm(N^2) 



% random samples removed
omega = sort(rPerm(1:nSamples));

% % %specific columns removed
% % omega = sort(rPerm);
% % omega = omega(1:nSamples);



Y = NaN(N)
Y(omega) = X(omega)


observations = X(omega)    % the observed entries
mu           = 1;        % smoothing parameter

% The solver runs in seconds
tic
Xk = solver_sNuclearBP( {N,N,omega}, observations, mu)
toc

sprintf('%10.1f',Xk')

% %     %Frobenius norm/L^2-norm
% %     fprintf('Relative error, no rounding: %.8f%%\n', norm(X-Xk,'fro')/norm(X,'fro')*100 );
% %     H (i)= norm(X-Xk,'fro')/norm(X,'fro')*100
% %     NoC (i) = nSamples %./ 10 %'/10' is for Columns
    
    %Manhattan norm/L1-norm
    fprintf('Relative error, no rounding: %.8f%%\n',norm(X-Xk,1)/norm(X,1)*100 );
    H (i)= norm(X-Xk,1)/norm(X,1)*100
    NoC (i) = nSamples %./ 10 %'/10' is for Columns


end
 H1(:,:,j) = H;
 NoC1(:,:,j)  = NoC;

end
% Converting the number stored in H1 into a string
H2 = num2str(H1)
NoC2 = num2str(NoC1)

% and change it into number again
Estimation_Error = str2num(H2)
Number_of_Samples = str2num(NoC2)



boxplot(Estimation_Error,Number_of_Samples)
xlabel('Number of Observable Samples');
ylabel('Estimation Error (%)');
title('Error between estimate and known solution');
% xlim([0 10]);
% ylim([0 100]);
grid;
% set(gca,'XTickLabel',[1:1:100])