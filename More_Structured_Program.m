clear;
clc;

%% More about structured programming
% 1. for loops
% 2. if-elseif-else conditonal statements
% 3. while loops
% 4. function calls

%% Creating functions in MATLAB
%clear;  % MatLab clears any variables in memory
%f = @(x) x^2;
%y = f(3)    % By leaving out the semi-colon this will be displayed


%% Algorithm 1: A recursive function for the factorial of n

% n = input('Number to get factorial of = ');
% 
% function fac = factorial(n)
%     if n == 0
%         fac = 0;
%     else
%         fac = n * factorial(n);
%     end
% end
% 
% result = fac(n)

%% Algorithm 2: Computation of e^x with a for-loop

% Error checking to make sure n is in range
% n =  input('Please input a number greater than or equal to 2: ');
% while n < 2
%     n =  input('Please input a number greater than or equal to 2: ');
% end
% 
% term = 1.0;
% y = term;
% 
% for k = 2:n-1
%     term = term*(1/(k-1));
%     y = y + term;
% end
% 
% result = y

%% Algorithm 3: Computation of e^x with a while loop
% 
% N = input('Please input a large natural number: ');
% x = input('Please input a real number: ');
% sigma = input('Please input a real number greater than 0: ');
% 
% % Error checking loops
% 
% while N < 100
%     N = input('Please input a LARGE natural number: ');
% end
% 
% while x < 0
%     x = input('Please input a REAL number: ');
% end
% 
% while sigma < 0
%     sigma = input('Please input a real number GREATER than 0: ');
% end
% 
% term = 1.0;
% y = term;
% diff = 2*sigma;
% k = 2;
% 
% while (diff > sigma && k < N)
%     term = (term * x) / (k - 1);
%     y = y + term;
%     diff = abs(term);
%     k = k + 1;
% end
% 
% result = y

%% Algorithm 4: Compute prime numbers less than a given natural number

% N = input("Please enter a natural number greater than 2: ");
% 
% while (N < 2)
%     N = input("Please enter a natural number GREATER than 2: ");
% end
% 
% primes = [2];
% length_primes = 1;
% 
% for n = 3:2:N
%     flag = 0;
%     for i = 1:length_primes
%         if mod(n, primes(i)) == 0
%             flag = 1;
%             break
%         end
%     end
%     if flag == 0
%         length_primes = length_primes + 1;
%         primes = [primes, n];
%     end
% end
% 
% result = primes



%% Algorithm 5: Demonstration of the for-else loop, and interating a list

% NEED TO MODIFY PYTHON TO MATLAB

% N = input("Please enter a natural number greater than 2: ");
% 
% while (N < 2)
%     N = input("Please enter a natural number GREATER than 2: ");
% end
% 
% primes = [2];
% length_primes = 1;
% 
% for n = 3:5:(N-1)
%     for x
% end

%% Algorithm 6: Recursive Fibonacci
% 
% N = input("Please enter a number greater than or equal to 0: ");
% y = 0;
% 
% while (N < 0)
%     N = input("Please enter a number GREATER THAN OR EQUAL TO to 0: ");
% end
% 
% function y = recursive_fibonacci(N)
%     if N == 1
%         y = 0;
%     elseif N == 2
%         y = 1;
%     else
%         y = recursive_fibonacci(N - 1) + recursive_fibonacci(N - 2);
%     end
% end
% 
% result = y

%% Algorithm 7: Improved Recursive Fibonacci

% function [y, z] = improved_recursive_fibonacci(N)
%     if N == 1
%         y = 0; z = 0;
%     elseif N == 2
%         y = 1; z = 0;
%     else
%         [a, b] = improved_recursive_fibonacci(N-1);
%         y = a + b;
%         z = a;
%     end
% end
%% Algorithm 8: An interesting recursion

% function L = Sierpinski(N, K)
%     L0 = [0 0;
%           1 0;
%           0.5 sqrt(3)/2];
% 
%     m0 = size(L0, 1);
% 
%     if N == 1
%         L = L0;
%     else
%         L = Sierpinski(N-1, K);
%         m = size(L, 1);
% 
%         temp = [];                  % Initialize tcc to be an empty array
% 
%         for i = 1:m0
%             x = ones(m,1) * L0(i,:); % Replicate vertex
%             temp = [temp; (L - x)/K + x];
%         end
%         L = temp;
%     end
% end

%% Algorithm 9

%% Algorithm 10
% function L = recursive_det(A)
% 
%     eps = 1e-10;
%     m = size(A, 1);
% 
%     if m == 1
%         y = A(1,1);
%         return;
%     end
% 
%     if abs(A(1,1)) > eps
%         k = A(1,1);
%         b = A(2:m, 1); % column below pivot
%         c = A(1,2:m)/k; % row excluding pivot
%         y = k * recursive_det(A(2:m,2:m) - b*c);
%     else
%         flag = 0;
%         for row = 2:m
%             if abs(A(row,1)) > eps
%                 A([1 row], :) = A([row 1],:); % Swap rows
%                 flag = 1;
%                 break;
%             end
%         end
% 
%         if flag == 1
%             k = A(1,1);
%             b = A(2:m,1);
%             c = A(1,2:m)/k;
%             y = -k * recursive_det(A(2:m,2:m) - b*c);
%         else
%             y = 0;
%         end
%     end
% end

%% Exercises

%% 1. Write a MATLAB program to compute n! using Algorithm 1.

% n = input('Number to get factorial of = ');
% 
% function fac = factorial(n)
%     if n == 0
%         fac = 0;
%     else
%         fac = n * factorial(n);
%     end
% end
% 
% result = fac(n)

%% 2. Consider the series {check notes} for -1 < x < 1.

% a. Write a pseudocode using the for-loop.

%   If we are using a for loop we must set a REAL end point NOT infinity
%   Because of this I will be using N as the end point for the for loop

%   Sequence is:    x + x^2/2 + x^3/3 + ... + x^N/N

%   N = user define
%   term = 0            % initialize to zero
%   y = term
%
%   for k = 0 to N do
%   term = x^k+1 / k+1
%   y = y + term


% b. Write a MATLAB program for the pseudocode in the previous item.
% Compare the result with the analytical sum -ln(1-x) for N = 10, N = 50,
% and N = 100.
% 
% x = input('Please input a number between -1 and 1: ');
% 
% while (x < -1 || x > 1)
%     x = input('Out of range, please input a number between -1 and 1: ');
% end
% 
% N = 10;                                 % Change this to change iterations
% 
% term = 0;                               % Initialize
% y = term;                               % =========
% 
% % Iterate through for loop
% for k = 0:N
%     term = (x^(k+1)) / (k+1);           % Update term
%     y = y + term;                       % Add term to previous term
% end
% 
% % Printing stuff
% result = y;                             % This is the value I just calculated
% ansum = -log(1 - x);                    % This is the true value I'm comparing to
% 
% disp('Analytical sum -ln(1-x):')
% disp(ansum)
% 
% disp('For loop implementation')
% disp(result)

% c. Re-write the pseudocode so that it is efficient.

%   N = user define
%   term = 0            % initialize to zero
%   y = term
%   a = k+1

%   for k = 0 to N do
%   term = x^a / a
%   y = y + term

% d. Rewrite the pseudocode using the while-loop.

%   For a while loop the end point is not defined, it can go forever unlike
%   a for loop. However, it would be better to set a designated top value
%   so that MATLAB does not... crash.

%   Sequence is:    x + x^2/2 + x^3/3 + ... + x^N/N

%   N = user define
%   term = 0            % initialize to zero
%   y = term
%   diff = 2sigma
%   k = 2
%
%   while  do
%   term = x^k+1 / k+1
%   y = y + term

% e. Write a MATLAB program for the pseudo-code in the previous item.
% Compare the result with the analytical sum -ln(1-x). You may use        
% sigma = 1e-5

%   Input for x with error checking
% x = input('Please input a number between -1 and 1: ');
% 
% while (x < -1 || x > 1)
%     x = input('Out of range, please input a number between -1 and 1: ');
% end
% 
% N = 50;                                 % Change this to change iterations
% sigma = 1e-5;                           % From Canvas
% 
% term = 0;                               % Initialize
% y = term;                               % =========
% 
% diff = 2*sigma;
% k = 0;
% 
% a = k + 1;
% 
% while (diff > sigma && k < N)
%     term = x^a / a;                     % Update term
%     y = y + term;                       % Add term to previous term
%     diff = abs(term);
%     a = a + 1;                          % Increment a
% end
% 
% % Printing stuff
% result = y;                             % This is the value I just calculated
% ansum = -log(1 - x);                    % This is the true value I'm comparing to
% 
% disp('Analytical sum -ln(1-x):')
% disp(ansum)
% 
% disp('While loop implementation')
% disp(result)

%% 3. Errors/Inefficient coding methods in mathematical computing

% a. Error in pseudocode: Does the pseudo-code presented in Algorithm 9
% below compute the first N terms of the series. If it is not correct how
% may you fix it?

% Algorithm 9 - Algorithm for sin(x)

% N = input("Please enter a natural value greater than 2: ");
% x = input("Please enter a real value for x: ");
% 
% while (N < 2)
%     N = input("Please enter a natural value GREATER than 2: ");
% end
% 
% sum = 1.0;
% term1 = x;
% term2 = 1.0;
% 
% for k = 1:N-1
%     term1 = -term1*x^2;
%     term2 = term2*(2*k)*(2*k+1);
%     sum = sum + term1/term2;
% end
% 
% result = sum

% b. Check the compuation of term1 in Algorithm 9. Is there a redundancy
% that you can remove?

% Yes, instead of saying term1 = x we can just put x in the loop.
% 
% N = input("Please enter a natural value greater than 2: ");
% x = input("Please enter a real value for x: ");
% 
% while (N < 2)
%     N = input("Please enter a natural value GREATER than 2: ");
% end
% 
% sum = 1.0;
% term2 = 1.0;
% 
% for k = 1:N-1
%     x = -x*x^2;
%     term2 = term2*(2*k)*(2*k+1);
%     sum = sum + x/term2;
% end
% 
% result = sum

% c. Write code for Algorithm 9 in MATLAB. With N = 15 and x = pi, step
% through the code in debug mode and observe the values of each of the
% variables. Are any variables getting extremely large? Is there a way you
% can rewrite the pseudo-code and code such that such large terms are
% avoided?

% N = 15;
% x = pi;
% 
% sum = 1.0;
% term2 = 1.0;
% 
% for k = 1:N-1
%     if isinf(sum)
%         error('Loop stopped at N = %d value cannot excede infinity', k)
%     end
%     x = -x*x^2;
%     term2 = term2*(2*k)*(2*k+1);
%     sum = sum + x/term2;
% end
% 
% result = sum

% d.

%% 4. 
% Compute the sum: 4[1 - 1/3 + 1/5 - 1/7 + ...] and compare the result with
% the analytical sum pi. Follow the steps of Problem 2. The above series is
% Leibniz formula for pi. Plot the partial sums to see the rate of
% convergence to pi. Now check convergence of the particular sums on the
% Machins formula for pi.

% Write the pseudocode using the for loop.

% Leibinz
% Input: N 
% sum = 0
% for k = 0 to N-1 do
%   term = ((-1)^k / (2*k + 1)
%   sum = sum + term
%
% result = 4 * sum
% Show result

% Machin
% Input: N
% sum1 = 0
% sum2 = 0
% for k = 0 to N-1 do
%   term1 = ((-1)^k * (1/5) ^ (2*k+1) / (2*k+1)
%   term2 = ((-1)^k * (1/239) ^ (2*k+1) / (2*k+1)
%   sum1 = sum1 + term1
%   sum2 = sum2 + term2
%
% result = 16 * sum1 - 4 * sum2
% Show result


% b. Write a MATLAB program for the pseudocode.

% Leibniz
% N = 10000;  % large number of iterations
% psums = zeros(1, N);
% 
% sumval = 0;         % Initialize sum, but we don't wanna use that name
% for k = 0:N-1
%     term = ((-1)^k) / (2*k + 1);
%     sumval = sumval + term;
%     psums(k+1) = 4 * sumval; % pi approx
% end
% 
% % Machins
% % Initialize a bunch of variables
% x1 = 1/5;
% x2 = 1/239;
% 
% atan1 = zeros(1, N);
% atan2 = zeros(1, N);
% 
% sum1 = 0;
% sum2 = 0;
% 
% for k = 0:N-1
%     term1 = ((-1)^k) *  x1^(2*k+1) / (2*k+1);
%     term2 = ((-1)^k) * x2^(2*k+1) / (2*k+1);
%     sum1 = sum1 + term1;
%     sum2 = sum2 + term2;
%     atan1(k+1) = sum1;
%     atan2(k+1) = sum2;
% end
% 
% machin = 16*atan1 - 4*atan2;
% 
% % Print stuff
% 
% fprintf('Leibniz seires: %.15f\n', psums(end))
% fprintf('Machin series: %.15f\n', machin(end));
% fprintf('Actual pi: %.15f\n', pi);
% 
% % Plot stuff
% M = 50;  % zoom window
% figure;
% hold on; grid on;
% 
% plot(1:M, psums(1:M), 'b', 'DisplayName', 'Leibniz series');
% plot(1:M, machin(1:M), 'k', 'LineWidth', 1.2, 'DisplayName', 'Machin series');
% yline(pi, 'r--', 'LineWidth', 1.5, 'DisplayName', 'π');
% 
% xlabel('Number of terms');
% ylabel('Partial Sum');
% title('Zoomed Convergence (First 50 Terms)');
% legend('show', 'Location', 'best');



%% 5. Error checking recursive Fibonacci.



%% 6. Use the tic and toc command to compare times for Algorithm 6
% N = 20;
% 
% tic();                              % Start timer
% y = recursive_fibonacci(N);         % Call algorithm 6
% t1 = toc();                         % End timer
% 
% tic();                              % Start timer
% [z, w] = improved_recursive_fibonacci(N);
% t2 = toc();
% 
% % Print stuff
% fprintf('Algorithm 6 time: %d\n', t1);
% fprintf('Algorithm 7 time: %d\n', t2);

%% 7. Modify algorithm 8

% function L = Sierpinski(N, K)
%     L0 = [cos(2*pi/5) -sin(2*pi/5);
%           1 0;
%           cos(2*pi/5) sin(2*pi/5);
%           -cos(pi/5) sin(4*pi/5);
%           -cos(pi/5) -sin(4*pi/5)];
% 
%     m0 = size(L0, 1);
% 
%     if N == 1
%         L = L0;
%     else
%         L = Sierpinski(N-1, K);
%         m = size(L, 1);
% 
%         temp = [];                  % Initialize tcc to be an empty array
% 
%         for i = 1:m0
%             x = ones(m,1) * L0(i,:); % Replicate vertex
%             temp = [temp; (L - x)/K + x];
%         end
%         L = temp;
%     end
% end
% 
% % Information given in question
% N = 5;
% K = 2.7;
% 
% L = Sierpinski(N, K);
% 
% % Plot
% figure;
% plot(L(:,1), L(:,2), '.k', 'MarkerSize', 4);
% axis equal;
% title(sprintf('Sierpinski recursion with N=%d, K=%d', N, K));
% xlabel('x');
% ylabel('y');

% Snowflakes :D

%% 8. Write pseudocode for computing twin primes less than a given #

% Pseudocode
% Input: N
% 
% // Find all primes less than N
% Initialize an empty list to store primes, lets call it P
% for i = 2 to N-1 do
%     isPrime = true
%     for j = 2 to floor(sqrt(i)) do
%         if i mod j == 0 then
%             isprime = false
%             break
%         end
%     end
%     if isPrime then
%         Append i to P
%     end
% end
% 
% // Find twin primes
% Initialize an empty list to store twin prime pairs, lets call it T
% 
% for k = 1 to length(P) - 1 do
%     if P[k+1] - P[k] == 2 then
%         Append (P[k], P[k+1]) to T
%     end
% end
% 
% Output T

% Implementation

% function twin_primes(N)
%     P = [];             % Initialize empty array for primes
%     for i = 2:N-1
%         isPrime = true;
%         for j = 2:floor(sqrt(i))
%             if mod(i,j) == 0
%                 isPrime = false;
%                 break;
%             end
%         end
%         if isPrime
%             P(end+1) = i;
%         end
%     end
% 
%     T = [];             % Initialize empty array for twin primes
%     for k = 1:length(P) - 1
%         if P(k+1) - P(k) == 2
%             T = [T; P(k), P(k+1)];
%         end
%     end
% 
%     disp('Twin primes less than N: ');
%     disp(T);
% end
% 
% twin_primes(100);

%% 9.

% see Algorithm 10
%% 10.

% a.
function c = poly_mult(a, b)
    n = length(a) - 1;
    m = length(b) - 1;
    r = n + m;
    c = zeros(1, r+1);

    for k = 0:r
        sumk = 0;
        for i = 0:k
            ai = 0; bi = 0;
            if i <= n
                ai = a(i+1);
            end
            if (k-i) <= m
                bi = b(k-1+1);
            end
            sumk = sumk + ai*bi;
        end
        c(k+1) = sumk;
    end
end

a = [-3 2 pi];
b = [1 0 3];
c = poly_mult(a, b);
disp(c);