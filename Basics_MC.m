clear;
clc;

%% Algorithm 1: While Loop Structure

% Initialize variables
% count = 1;                
% sigma = 1;
% 
% while sigma > 10^-5 & count <= N  % While these conditions are met
% 
%     count = count + 1             % Count increments by 1
%     sigma = sigma + 1             % Not sure what this was supposed to
%                                   % update to
% end

%% Algorithm 2: Fibonacci Sequence 1

% Initialize starting points
% a(1) = 0;               
% a(2) = 1;
% 
% for k = 3:N                       % From third position to Nth position              
%     a(k) = a(k-2) + a(k-1)        % Perform Fibonacci formula
% end

%% Algorithm 3: Fibonacci Sequence 2

% Not very necessary series of loops
% if N == 1
%     a(1) = 0;
% elseif N ==2
%     a(1) = 0;
%     a(2) = 1;
% else
%     a(1) = 0;
%     a(2) = 1;
% end
% 
% for k = 3:N                       % From third position to Nth position              
%     a(k) = a(k-2) + a(k-1)        % Perform Fibonacci formula
% end


%% Algorithm 4: Unknown Program

% This code is computing a factorial? N is predefined within my code to be
% 3 and the factorial of 3 is 6
% x = 0;
% for k = 1:N
%     x = x + k
% end

%% Algorithm 5: Recaman Sequence
% 
% Error checking to make sure N is in range
% N = input('Please insert a natural number greater than or equal to 3: ');
% while N < 3
%     N = input('Please insert a natural number greater than or equal to 3: ');
% end
% 
% Following the Pseudocode for the flags
% a(1) = 0;
% for k = 2:N
%     temp = a(k-1) - (k-1);
%     flag = 0;
% 
%     if temp < 0
%         flag = 0;
%     else
%         if isempty(temp)
%             flag = 1;
%         end
%     end
% 
%     if flag == 0
%         a(k) = a(k-1) + (k-1)
%     else
%         a(k) = temp
%     end
% end

%% Exercises

%% 1. Create the matrix:
%       [2 -3 -1 4]
%   A = [0 -1  1 2]

% A = [2, -3, -1, 4; 0, -1, 1, 2]
% 
% % a. Add column Ab between second and third columns of A, call the result B
% Ab = [0; 0];
% 
% B = [A(:,1:2) Ab A(:,3:end)]
% 
% % b. Delete the second column of B and call the result C
% B(:, 2) = [];
% C = B


%% Table 1 from Basics.pdf

% A. The Federal Reserve has a target annual inflation rate of 2%, and it
% uses monetary policy to keep infation in check and stablize the economy
% when inflation rises above that benchmark. Compute the average inflation
% rate for the period 2001 - 2022.

IF = [8.9, 3.8, 3.8; 3.9, 3.8, 1.1; 4.4, 4.4, 4.6; 
      6.1, 3.1, 2.9; 2.7, 2.7, 2.5; 3.3, 1.7, 1.6; 
      2.7, 3.4, 1.6; 2.4, 1.9, 3.3; 3.4, 2.5, 4.1; 
      0.1, 2.7, 1.5; 3.0, 1.7, 1.5; 0.8, 0.7, 2.1;
      2.1, 1.9, 2.3; 1.4, 7.0, 6.5];

% Turn this into a array for calculation

r = reshape(IF.',1,[]);
%N = 41;

%for k = 1:N
%    x(:, k) = (1 + (r(k)/100));
%end

%rbar = 100*(exp((log(x)/N)-1));

Y = 2022 - 1981
Yn = 2002 - 1981
YEARS = 2022-2002;

%A = sum(rbar(:, YEARS:end),2)/YEARS

A = sum(r(:, 22:end),2)/21

% B. TThe last two were the COVID years. Compute the average inflation rate
% for the period 2001 - 2020

A = sum(r:, 22:end-2), 2)/19