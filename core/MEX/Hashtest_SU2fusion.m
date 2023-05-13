keylist = {};
for R1 = 1:80
    for R2 = 1:80
        for R = ((R1-R2)+1):(R1+R2-1)
            keylist{end+1} = char([R1,R2,R]);
        end
    end
end

tic;
bucketlist = Hashtest(keylist);
toc;

tic;
shitty_map = containers.Map('KeyType','char','ValueType','double');
for i = 1:length(keylist)
    shitty_map(keylist{i})=1;
end
toc;
