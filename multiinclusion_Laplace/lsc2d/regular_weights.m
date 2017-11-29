function [x,w] = regular_weights(a)

par = [2;3;4;5;7;10];
for i = 1:length(par)
    if par(i)==a
        key = i;
    end
end

lenth = [2;3;4;6;8;12];
xp = load('nodes_regular.dat');
if key==1
    starting = 1;
else
    starting = sum(lenth(1:key-1))+1; 
end

x = xp( starting:starting+lenth(key)-1,1);
w = xp(starting:starting+lenth(key)-1,2);