function out = SITE_sector_exist(site, secrep)
%SITE_SECTOR_EXIST Summary of this function goes here
%   Detailed explanation goes here
key = char(secrep);
out = site.multiplet_dims.isKey(key);

end

