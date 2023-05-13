function data = read_FCIDUMP_AB(infile)
%READ_FCIDUMP_AB Summary of this function goes here
%   Detailed explanation goes here
    fid = fopen(infile);
    elems = {0};
    next = false;
    while(~isequal(elems{end},'&END'))
        tline = fgetl(fid);
        elems = split(tline,{',',' '});
        for i = 1:length(elems)
            if next
                NORB = str2double(elems{i}); 
            end
            if isequal(elems{i},'NORB=')
                next = true;
            else
                next = false;
            end
        end
    end
    T = zeros(NORB,NORB);
    V = zeros(NORB,NORB,NORB,NORB);
    Vnum=0;
    Tnum=0;
    while ~feof(fid)
        tline = fgetl(fid);
        tmp = str2num(tline);
        if isequal(tmp(2:5),[0,0,0,0])
            Eshift = tmp(1);
        elseif isequal(tmp(4:5),[0,0])
            T(tmp(2),tmp(3)) = tmp(1);
            %Tnum = Tnum + 1;
        else
            V(tmp(2),tmp(3),tmp(4),tmp(5)) = tmp(1);
            %Vnum = Vnum + 1;
        end
    end
    %disp([Tnum, Vnum]);    
    data = struct();    
    data.norb = NORB;
    data.T = T;
    data.V = V;
    data.Eshift = Eshift;
    fclose(fid);

end

