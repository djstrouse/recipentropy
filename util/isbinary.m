function [flag] = isbinary(x)

% Checks if input is binary (i.e. elements are all zero or one)

if sum(((x(:)==0)+(x(:)==1))~=1)>0
  flag = false;
else
  flag = true;
end

end

