function expt=filtesweeps(expt,bOR,varargin)
[expt.sweeps ind]=filtsweeps(expt.sweeps,bOR, varargin{:});
if isfield(expt.wc,'data');
    expt.wc.data=expt.wc.data(ind,:);
end
if isfield(expt.wc,'lowpassdata')
    expt.wc.lowpassdata=expt.wc.lowpassdata(ind,:);
end
if isfield(expt.wc,'highpassdata')
    expt.wc.highpassdata=expt.wc.highpassdata(ind,:);
end
if isfield(expt.wc,'normdata')
    expt.wc.normdata=expt.wc.normdata(ind,:);
end
end
