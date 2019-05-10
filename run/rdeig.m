%     Read Hes%03i

clear
clc
close all

mksz=8;

fn = [999]; %1:10;

h1=figure(1);

for i=fn
  fname = sprintf('%s%3.3i','Eig',i);
  ain = importdata(fname);
  [r c] = size(ain);
  Ar = ain(:,1:2:c);
  Ai = ain(:,2:2:c);
  A = Ar + 1i*Ai;

  plot(Ar,Ai, '.', 'MarkerSize',mksz); hold on
  
end  


