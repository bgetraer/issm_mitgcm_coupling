function s = appNum(i,n)
% s = appNum(i,n)
% i = integer to be converted
% n = length of string (incl leading zeros)

if (round(i)~=i)
 error('input must be integer');
end

a = num2str(i);

if (length(a)>n)
 a=a(end-(n-1):end);
elseif (length(a)<n)
 m = n-length(a);
 aa = '0';
 for k=1:m-1
  aa=[aa '0'];
 end
 a = [aa a];
end

s=a;
return
