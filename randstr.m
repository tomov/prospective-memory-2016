function string = randstr( len )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

SET = char(['a':'z' '0':'9']) ;
NSET = length(SET) ;

N = len ; % pick N numbers
i = ceil(NSET*rand(1,N)) ; % with repeat
R = SET(i) ;
string = R;
end

