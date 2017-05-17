function int_list = str_to_int_list(str_list)
int_list = [];
last = false;
while ~last
    ic = strfind(str_list, ',');
    if isempty(ic)
        ic = length(str_list) + 1;
        last = true;
    end
    str_item = str_list(1 : ic - 1);
    str_list = str_list(ic + 1 : end);
    ic = strfind(str_item, '-');
    if isempty(ic)
        int_item = [str2num(str_item)];
    else
        strt = [str2num(str_item(1 : ic - 1))];
        stop = [str2num(str_item(ic + 1 : end))];
        int_item = strt : stop;
    end
    int_list = [int_list int_item];    
end
end