function out = SITE_get_site_irreps(obj)
            out = cellfun(@(x) double(x), obj.sector_multiplets.keys(),'UniformOutput',false);
end

