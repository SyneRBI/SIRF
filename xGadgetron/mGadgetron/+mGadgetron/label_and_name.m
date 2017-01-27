function [label, name] = label_and_name(input)
name = strtrim(input);
i = strfind(name, ':');
if isempty(i)
    label = '';
else
    if i == 1
        label = '';
    else
        label = strtrim(name(1 : i - 1));
    end
    name = strtrim(name(i + 1 : end));
end
end