function fileProcess(path,suffix)
files = dir([path,'*.',suffix]);
len = length(files);
for i = 1:len
    filePath = [path,files(i).name]; %得到文件路径
    %文件处理部分
end