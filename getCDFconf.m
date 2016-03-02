function confint = getCDFconf (data,percent)

% get upper and lower confidence bounds of non-normal data
% percent is percent of data within confidence bounds (input the upper
% bound)
confint=[]; %confint returns upper and lower confidence boud
% for example, if percent = 85, confint is lower 15% and upper 85%
upper= percent/100;
lower = (100-percent)/100;
[xb,pb]=empcdf(data);
tmp1=find(pb<=upper);
tmp2=find(pb>=lower);
tmp=intersect(tmp1,tmp2);
confint(1)=xb(min(tmp));
confint(2)=xb(max(tmp));
