% close all;
% clear;
% clc;
% temp = 0

N   = 10;       % the matrix is N x N
r   = 2;        % the rank of the matrix
% df  = 0:0.1:0.9;
df  = linspace(0,0.9,10)
Number_of_Experiments = 1


% for point_row = 1:1:Number_of_Experiments
% for point_column = 1:1:N   
for j = 1:Number_of_Experiments   
iMax    = 10;
X       = randi(iMax,N,r)*randi(iMax,r,N) % Our target matrix
rPerm   = randperm(N^2) % use "randsample" if you have the stats toolbox
    
    for i = 1:length(df)
        %nSamples    = 2*df; % number of observed entries
        nSamples = (N^2)-(df(i)*(N^2));

        % For this demo, we will use a matrix with integer entries
        % because it will make displaying the matrix easier.
        % iMax    = 5;
        % X       = randi(iMax,N,r)*randi(iMax,r,N) % Our target matrix
        % rPerm   = randperm(N^2) % use "randsample" if you have the stats toolbox

%         %random samples removed
%         omega = sort(rPerm(1:nSamples));
        
%         %random columns removed
%         k = randperm(N)
%         M = NaN(k)
%         omega = X(:,k(1:(nSamples/10)));
        
        %specific columns removed
        omega = sort(rPerm);
        omega = omega(1:nSamples);

        Y = NaN(N);
        Y(omega) = X(omega);

disp('The "NaN" entries represent unobserved values');
disp(Y)
        observations = X(omega);    % the observed entries
        mu           = .01;        % smoothing parameter

    % The solver runs in seconds
    tic
    Xk = solver_sNuclearBP( {N,N,omega}, observations, mu );
    toc
    
%     %Frobenius norm/L^2-norm
%     fprintf('Relative error, no rounding: %.8f%%\n', norm(X-Xk,'fro')/norm(X,'fro')*100 );
%     H (i)= norm(X-Xk,'fro')/norm(X,'fro')*100
%     NoC (i) = nSamples /100%./ 10 %'/10' is for Columns
    
    %Manhattan norm/L1-norm
    fprintf('Relative error, no rounding: %.8f%%\n',norm(X-Xk,1)/norm(X,1)*100 );
    H (i)= norm(X-Xk,1)/norm(X,1)*100
    NoC (i) = nSamples %./ 10 %'/10' is for Columns
   
    end
    
        H1(:,:,j) = H;
        NoC1(:,:,j)  = NoC;
end

H2 = num2str(H1)
NoC2 = num2str(NoC1)

Estimation_Error = str2num(H2)
Number_of_Samples = str2num(NoC2)

% temp = temp + 1
% Estimation_Error(temp) = H1(point_row, point_column)
% Number_of_Samples(temp) = NoC1(point_row, point_column)
% Estimation_Error((point_row-1)*N+point_column) = H1(point_row, point_column)
% Number_of_Samples((point_row-1)*N+point_column) = NoC1(point_row, point_column) 

% end 
%   
% end

% disp('Recovered matrix (rounding to nearest .0001):')
% disp( round(Xk*10000)/10000 )
% % and for reference, here is the original matrix:
% disp('Original matrix:')
% disp( X )

% The relative error (without the rounding) is quite low:
% fprintf('Relative error, no rounding: %.8f%%\n', norm(X-Xk,'fro')/norm(X,'fro')*100 );

% boxplot(Estimation_Error,Number_of_Samples)
% 
% 
% xlabel('Number of Observable Samples');
% ylabel('Estimation Error (%)');
% title('Error between estimate and known solution');
% % xlim([0 100]);
% % ylim([0 100]);
% grid on;
% % set(gca,'XTickLabel',[1:1:100])
% 
