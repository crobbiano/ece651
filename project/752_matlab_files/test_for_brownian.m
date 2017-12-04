rn=drg-g;
noise=u_x0_noisy-u_x0;
% do moving window derivitive and then plot by bins

% for i=1:length(noise)-1
%    temp(i) = (noise(i+1)-noise(i))/tstepsize;
% end
% for i=1:length(noise)-2
%    temp2(i) = (noise(i+2)-noise(i))/(2*tstepsize);
% end
% for i=1:length(noise)-3
%    temp3(i) = (noise(i+3)-noise(i))/(3*tstepsize);
% end
% figure (7)
% subplot(3,1,1);
% [h,pd]=histfit_w_out(temp); % mu=-0.198154987977091;sigma=4.832938893240116
% str=sprintf('Width=1; Mu=%f; Var=%f;',pd.mu, pd.sigma);
% title(str);
% subplot(3,1,2);
% [h,pd]=histfit_w_out(temp2) % mu=-0.139430362966122; sigma=3.136994536490363
% str=sprintf('Width=2; Mu=%f; Var=%f;',pd.mu, pd.sigma);
% title(str);
% subplot(3,1,3);
% [h,pd]=histfit_w_out(temp3) % mu=-0.079451195534108; sigma=1.637753594314555
% str=sprintf('Width=3; Mu=%f; Var=%f;',pd.mu, pd.sigma);
% title(str);

for i=1:length(rn)-1
   temp(i) = (rn(i+1)-rn(i))/tstepsize;
end
for i=1:length(rn)-2
   temp2(i) = (rn(i+2)-rn(i))/(2*tstepsize);
end
for i=1:length(rn)-3
   temp3(i) = (rn(i+3)-rn(i))/(3*tstepsize);
end
figure (6)
subplot(3,1,1);
[h,pd]=histfit_w_out(temp); % mu=-0.198154987977091;sigma=4.832938893240116
str=sprintf('Width=1; Mu=%f; Var=%f;',pd.mu, pd.sigma);
title(str);
subplot(3,1,2);
[h,pd]=histfit_w_out(temp2) % mu=-0.139430362966122; sigma=3.136994536490363
str=sprintf('Width=2; Mu=%f; Var=%f;',pd.mu, pd.sigma);
title(str);
subplot(3,1,3);
[h,pd]=histfit_w_out(temp3) % mu=-0.079451195534108; sigma=1.637753594314555
str=sprintf('Width=3; Mu=%f; Var=%f;',pd.mu, pd.sigma);
title(str);