files = dir('*txt');
numFiles = length(files);

for i = 1:numFiles
    d = load(files(i).name);
    data = d';
    save(['C' files(i).name(1:end-4) '.mat'],'data');
end