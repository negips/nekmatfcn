%     Read Hes%03i

clear
clc

fn = 1:142;

for i=fn
  fname = sprintf('%s%3.3i','Hes',i);
  ain = importdata(fname);
  [r c] = size(ain);
  Ar = ain(:,1:2:c);
  Ai = ain(:,2:2:c);
  A = Ar + 1i*Ai;
  
  fA = logm(A);
  resid(i) = abs(fA(r,1));
end  

[fn' resid']
semilogy(resid)

