function[] = construct_Table(K11,N11,error_max)
%determining convergence rate
con_rate = zeros(size(K11,2),size(N11,2));
for i = 1:size(K11,2)-1
  for j = 1:size(N11,2)
   con_rate(i+1,j) = log(error_max(i,j)/error_max(i+1,j))/log(2);
  end
end

disp('CONVERGENCE RATES')

N1 = error_max(:,1);
N2 = error_max(:,2);
N3 = error_max(:,3);
N4 = error_max(:,4);

Rate1 = con_rate(:,1);
Rate2 = con_rate(:,2);
Rate3 = con_rate(:,3);
Rate4 = con_rate(:,4);

no_of_elements = {'K = 5 ','K  = 10 ','K = 20','K = 40','K = 80'};
T = table(N1,Rate1,N2,Rate2,N3,Rate3,N4,Rate4,'Rownames', no_of_elements)