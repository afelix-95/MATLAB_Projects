function decoderdrive
%Alison Felix de Araujo Maia, CAAM 210, Fall 2015, Lab X
%decoderdrive.m
%Description: This function decode encoded messages using the Metropolis
%Algorithm.
for i = 1:3 %The for-loop is to make sure that all 3 messages appear well decoded at least once.
text1 = decoder('encodedtext1.txt');
text2 = decoder('encodedtext2.txt');
text3 = decoder('encodedtext3.txt');
disp('Text 1:');
disp(text1);
disp(' ');
disp('Text 2:');
disp(text2);
disp(' ');
disp('Text 3:');
disp(text3);
disp(' ');
end
return

function text = decoder(filename)
T = fileread(filename);
T = double(T);
for j = 1:length(T)%for-loop to change the ascii index numbers in numbers from 1 to 27 (a-z + space).
    T(j) = downlow(T(j));
end
y = randperm(27);%Random initial guess.
for j = 1:10000
    k1 = randi(27);
    k2 = randi(27);
    ymaybe = y;
    ymaybe([k1,k2]) = ymaybe([k2,k1]); %Here is generated the ymaybe by switching two random elements picked with the ''randi'' function above.
    l1 = loglike(T, y);
    l2 = loglike(T, ymaybe);
    if l1 < l2 %Conditions to replace the first guess (y) by the new guess (ymaybe) based on the log-liklihood calculated above.
        y = ymaybe;
    elseif rand < exp(-l1+l2)
        y = ymaybe;
    end
end
for j = 1:length(T) %for-loop to change the numbers inside the array (1 through 27) back to ascii index numbers.
    T(j) = y(T(j));
    T(j) = downlowinv(T(j));
end
text = char(T);%change the vector from ascii index numbers back to characters. 
return

function x = downlow(n)%Function used to turn the ascii index numbers into numbers from 1 through 27
if n == 96
    x = 27;
else
    x = n-96;
end
return

function c = downlowinv(n)%Function used to turn the numbers from 1 through 27 into ascii index numbers 
if n == 27
    c = 32;
else
    c = n+96;
end
return

function x = loglike(T, y)%Function used to calculate the log-liklihood
x = 0; M = textread('letterprob.mat');
for j = 1:(length(T)-1)
    x = x + log(M(y(T(j)),y(T(j+1))));
end
return



