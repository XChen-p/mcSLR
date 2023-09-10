function out = thresh(input, c, mask)

if nargin == 2
    mask    =   input;
end

out =   input.*(abs(mask)>c);
