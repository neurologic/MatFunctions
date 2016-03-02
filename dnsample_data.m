function data=dnsample_data(data,oldSampRate,newSampRate)
%for this to work, the data must be along one row, with columns as samples
% didflip=0;
% if find(max(size(data)))==1
%     data=data';
%     didflip=1;
% end

bin=round(oldSampRate/newSampRate);

data=data(:,1:bin:end,1);
% if didflip==1
%     
%     data=data';
% end
