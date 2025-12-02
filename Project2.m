% Thomas Haimes
% Math 4330
% Project 2 - Minimum Norm Least Squares

%Clear console for neatness
clc;
clear;
close all;

%% Inputs

% A and b for manual input
fprintf("PLEASE FOLLOW THE FOLLOWING INSTRUCTIONS FOR MANUAL INPUT");
A = input('\nEnter the coefficient matrix of A (as [a1, a2, a3; b1, b2, b3; etc] format): ');
b = input('\nEnter the coefficient matrix of b: ');

% Test case for QR
% A = [1,2,0; 2,1,2; 1,1,-3];

% A and b for exercises, uncomment the section and comment out above ^
% 1.
% A = [4,0,0,0; -4,8,-4,2; -7, 7, -3, 3; -2, 2, -2, 4];
% b = [1;2;3;4];

% 2.
% A = [1,2;1,4;1,9];
% b = [1;2;3];

% 3.
%  A = [1, 12800;     % Values of x1 + x2*v
%       1, 13100;
%       1, 12950;
%       1, 13580;
%       1, 13760;
%       1, 14100;
%       1, 13890;
%       1, 14320;
%       1, 14880;
%       1, 15220];
% 
% 
%  xbar = [; ]; % Empty array to store new values for xbar
% 
%  b = [857000; 917000; 841750; 869120; 894400;       % Values of y
%       874200; 972300; 887840; 967200; 942400];
% 
%  x =LSQR(A, b)
%  v = linspace(1.2e4,1.6e4,51);
%  y = x(1) +x(2)*v;
% plot(v,y,'b');
% hold on;
% grid on;
% plot(A(:,2),b,'r*')
% LSQRN(A,b )
 %simpleQR(A);

 % A = [0, 0; 0 , 0];
 % b = [1; 2];
 % x = LSQR(A, b)

%% Algorithm 1 - Simple QR decomposition algorithm
function [Q, R] = simpleQR(A)
  %First initializes Q and R to the m x n and n x n zero matrices
  [m, n] = size(A);
  R = zeros(n, n);
  Q = zeros(m, n);
 
  % Compute the first column of Q and R
  R(1,1) = norm(A(:, 1));       % Magnitude of the first column of A
  Q(:,1) = A(:, 1) / R(1,1);    % Normalize to form first Q column
 
  % Iterate for remaining columns
  for k = 2:n
    vector_proj = zeros(m, 1); % Initialize projection vector
    for i = 1:(k-1)
        R(i, k) = Q(:,i)'*A(:, k); % Compute projection coefficient
        vector_proj = vector_proj + R(i, k)*Q(:, i); % Accumulate projections
    end
    residual = A(:, k) - vector_proj; % Orthogonal component
    R(k,k) = norm(residual);          % Normalize residual to get diagonal element
    Q(:,k) = residual / R(k,k);       % Normalize to form next column of Q
  end

  % Display stuff
  fprintf("QR decomposition of A:\n");
  disp("QR:");
  disp([A Q*R]);
  disp("Q'Q:");
  disp(Q'*Q);
  disp("Q:");
  disp(Q);
  disp("R:");
  disp(R); 
end

 %simpleQR(A)

 %% Algorithm 2 - Computation of the solution to the least squares problem
 
 function x = LSQR(A, b)
    [Q, R] = simpleQR(A);       % Compute Q and R
    b = Q'*b;                   % Transform b into Q-space
    x = backward(R, b);         % Solve upper-triangular system

    fprintf("Least squares of A: \n");
    disp(x);
 end

%LSQR(A, b)

 %% Algorithm 3 - Backward substitution algorithm

 function x = backward(U, z)
    n = length(U);                          % Determine size of system
    x = zeros(n, 1);                        % Initialize x
    x(n) = z(n) / U(n, n);                  % Solve last equation directly
    for i = n-1:-1:1                        % Loop backwards through rows
        sum = 0;
        for j = i+1:n
            sum = sum + (U(i, j) * x(j));   % Accumulate known terms
        end
        x(i) = (z(i) - sum) / U(i, i);      % Solve for current x(i)
    end

 end


 %% Algorithm 4 - Improved QR decomposition algorithm when columns of A are not LI
 function [Q, R, rank] = QR(A)
    esp = 1e-12; % Small number for tolerance
    [m,n] = size(A);
    % Initialize Q and R to be mxn and nxn zero matrices
    R = zeros(n,n);
    Q = zeros(m,n);
    rank = 0; 

    for j = 1:n
        vector_proj = zeros(m, 1);  % Initialize projection for column j
        for i = 1:rank
            vector_proj = vector_proj + R(i,j)*Q(:,i); % Subtract existing projections
        end

        residual = A(:,j) - vector_proj;    % Othogonal component
        norm_residual = norm(residual);     % Compute nrom of residual

        % Only add a new Q column if residual is significant
        if norm_residual > esp
            rank = rank + 1;                        % Increment rank
            R(rank, j) = norm_residual;             % Store new R diagonal
            Q(:, rank) = residual/norm_residual;    % Normalize new Q column

            % Only add a new Q column if residual is significant
            for k = j+1:n 
                R(rank, k) = Q(:, rank)' * A(:,k);
            end
        end
    end

    % Truncate Q and R to computed rank
    Q = Q(: ,1:rank);
    R = R(1:rank, :);

  % Display stuff
  fprintf("Improved QR decomposition of A:\n");
  disp("QR:");
  disp([A Q*R]);
  disp("Q'Q:");
  disp(Q'*Q);
  disp("Q:");
  disp(Q);
  disp("R:");
  disp(R);
 end

 %QR(A)

 %% Algorithm 5 - Forward substitution algorithm
 function z = forward(L, b)
    n = length(b);  % Initialize n to the length of the vector so it exists
    z = zeros(n,1); % Initialize z to be the zero nx1 vector
   
    z(1) = b(1)/L(1,1); % Solve for first element
    for i = 2:n         % Loop from 2 to n
        sum = 0;        % Initialize sum
        for j = 1:i-1   % Compute dot product of row with known z's
            sum = sum + L(i,j)*z(j);
        end
        z(i) = (b(i)-sum)/L(i,i); % Solve for z(i)
    end
 end

 %% Algorithm 6 - Minimum norm solution to least squares problem

 function x = LSQRN(A, b)
    [Q,R,rank] = QR(A);
    b = Q'*b; % Find the number of columns of the A matrix and call it n
    [m, n] = size(A);
    if rank == n
    % If rank is n, then there is no need to find the minimum norm solution
    % There is only one unique solution in this case
        x = backward(R, b);
    else
    % In this case, rank < n
    % The least squares solution is not unique
    % We have to do another QR decomposition to find the minimum-norm least
    % squares solution
        [Q1, R1, rank1] = QR(R');
        z = forward(R1', b);
        x = Q1*z;
    end

    fprintf("Least squares of A: \n");
    disp(x);
 end

%% Linear Regression
% From Canvas
% mat = readmatrix('test.csv');   % Read in data set
% I = randi([1, 6700], 1, 6000);      % Randomly select upt to 6000 homes
% J = unique(I);                      % Modeling set
% K = setdiff(1:6700, J);             % Test set
% A = mat(J, 3:end);                  % Create the A matrix
% b = mat(J, 2);                      % Create the b vector
% x = lsqr(A, b);                     % Create a linear model
% 
% % Test
% [mat(K,3:end)*x mat(K,2)]

%% Second set of exercises

% Test case
% A = [1,5,9; 2,6,10; 3,7,11; 4,8,12]
% b = [1;2;3;4]

%A = [1, 2, 3; 2, 4, 6];
%b = [1; 1];

% A = [0, 1, 5, 9, -1; 0, 2, 6, 10, 1; 0, 3, 7, 11, 1; 0, 4, 8, 12, 1];
% b = [1; 2; 3; 4];

LSQR(A, b); % Works
LSQRN(A, b); % Doesn't work keep debugging