function load_metatoes

metatoes = load('metatoes.mat');
metatoesVars = metatoes;
vars = fieldnames(metatoes);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= metatoesVars.' vars{ifield} ';'];
    eval(s)
end
