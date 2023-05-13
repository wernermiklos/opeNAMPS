function out = SITE_get_tau_sectors(obj)
        out = cell(1,length(obj.sector_multiplets));
        i = 1;
        for keycell = obj.sector_multiplets.keys
            key = keycell{1};
            irrep = double(key);
            sectordim = obj.sector_multiplets(key);
            out{i} = {{irrep},sectordim};
            i = i + 1;
        end
end

