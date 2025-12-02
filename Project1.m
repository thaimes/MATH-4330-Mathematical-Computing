% Thomas Haimes
% MATH 4330
% Project 1 - RSA Algorithm

%Clear console for neatness
clc;
clear;
close all;

%% Basics of modern cryptography

%% Algorithm 1: Extended Euclidean algorithm for gcd of a and p

% Testing values
% a = 63;  % Numerator
% p = 11;  % Denominator

function [oldr, olds] = gcd_alg(a, p)
    r = [a;p];                      % Initialize remainder
    s = [1;0];                      % Initialize s / Bezout coeff a
    t = [0;1];                      % Initialize t / Bezout coeff b

    while (r ~= 0)                  % While remainder is not 0
        q = floor(r(1)/r(2));       % Round result to nearest int
        r = [0, 1; 1, -q] * r;      % Update remainder using matrix mult
        s = [0, 1; 1, -q] * s;      % Quotient by gcd
        t = [0, 1; 1, -q] * t;      % Quotient by gcd
    end

    if s(1) < 0
        s(1) = s(1) + p;
    end

    oldr = r(1);
    olds = s(1);
end
% 
% gcd(a, p)

%% Algorithm 2: Algorithm for multiplicative inverse modulo p
% Test values
% a = 3;
% p = 6;

function [s, flag] = inverse_mod_n(a, p)
    [r,s] = gcd_alg(a, p);          % Extended Euclidean Algorithm
    if r ~= 1                       % If remainder not 1
        flag = false;               % False flag
    else
        flag = true;                % True flag
    end
end

%inverse_mod_n(a, p)


%% Exercises for Modular Arithmetic

%%1. 
% a = 0:20;
% p = 21;
% 
% 
% for k = 1:p-1
%     gcd(a(k), p);
% end

%%2. On OneNote

%% RSA Algorithms

%% Algorithm 1 - Exponential Modulo N Algorithm

function r = ModExp(m, e, N)
    a = dec2bin(e)-'0'; % Convert e to binary
    r = 1;              % Required result
    t = m;              % Powers of m

    n = length(a);

    for i = n:-1:1
        if a(i) == 1
            r = mod(r*t,N);
        end
        t = mod(t*t,N);
    end
end

% m = "message";
% e = 11;
% N = 15;

% ModExp(m, e, N)

%% Algorithm 2 - RSA algorithm simulation

function RSA_sim()
    % https://github.com/srmalins/primelists/blob/master/100primes.txt
    primeslist = load("primes.txt");                % Load list of primes
    length_primes = length(primeslist);             % Number of primes
    
    % Primes selected randomly from primeslist
    prime1 = primeslist(randi([2, length_primes])); % p randomly generated
    prime2 = primeslist(randi([2, length_primes])); % q randomly generated
    %prime1 = 3;                                     % p - exercise 2
    %prime2 = 7;                                     % q - exercise 2

    % prime1 and prime2 must be distict
    while prime1 == prime2
        prime2 = randi([2, length_primes]); % regenerate q if q = p
    end
    
    % Find the least common multiple (LCM) of (prime1 - 1)(prime2 - 1)
    [g, ~] = gcd_alg(prime1 - 1, prime2 - 1);  
    lambda = ((prime1 - 1) * (prime2 - 1)) / g;

    % Public Key (e)
    found = false;                  % flag to say if the key is found
    while ~found
        e = randi([3, lambda - 1]); % e - randomly generated between 3 and lambda
        %e = 5; % e - defined in exercise 2
        [k, found] = inverse_mod_n(e, lambda); % k is Alice's private key
    end

    N = prime1 * prime2;                    % Modulus N for public key

    % Print key information for debug
    fprintf('Public key (N, e) = (%d, %d)\n', N, e);
    fprintf('Private key (d) = %d\n', k);

    % Get message
    m = input('Enter message to encrypt here: ', 's'); % User enters string
    msglen = length(m);                                % Length for future calc

    % Encryption time
    encoded = zeros(1, msglen); % Populate whole array with 0s

    % For loop to step through entire message
    for i = 1:msglen
        x = double(m(i));       % Convert each character to numeric
        %x = double('H');
        encoded(i) = ModExp(x, e, N); % Perform modular exponentiation
    end

    % Print encrypted message values... for not reading
    fprintf('Encoded numeric message: %d\n', encoded);

    % Decryption time
    decoded = char(zeros(1, msglen)); % Populate whole array with 0s

    % For loop to step through entire message
    for i = 1:msglen
        x = ModExp(encoded(i), k, N);   % Perform modular exponentiation w/ private key
        decoded(i) = char(x);           % Convert each value to characters
    end

    % Print decrypted concatenated message for reading
    fprintf('Decoded message: %s\n', decoded);
end

% Run simulation
% Whatever message you want to encrypt will be typed into the console
RSA_sim