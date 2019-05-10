%     Read Hes%03i

clear
clc
close all

fn = [2:19]; %1:142;

for i=fn
  fname = sprintf('%s%3.3i','hessenberg',i);
  ain = importdata(fname);
  [r c] = size(ain);
  A = ain(1:c,:);

  e = eig(A);
  plot(real(e),imag(e), '.'); hold on
end  

%[fn' resid']
%semilogy(resid)

