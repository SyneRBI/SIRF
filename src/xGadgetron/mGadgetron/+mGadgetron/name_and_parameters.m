function [name, par] = name_and_parameters(input)
name = strtrim(input);
i = strfind(name, '(');
if isempty(i)
    par = [];
else
    j = strfind(name(i + 1 : end), ')');
    par = strtrim(name(i + 1 : i + j - 1));
    if i == 1
        name = '';
    else
        name = strtrim(name(1 : i - 1));
    end
end
end