function [v, u, a] = give_the_weights(q, ker)

switch ker
    case 0
    xp = load('nodesr.dat');
    lenth = [2;4;8];
    par =[2;4;7];
    case 1
    xp = load('nodes_sqrtx.dat');
    lenth = [4;8;16];
    par = [3;5;10];
    case 2
    xp = load('nodes_logx.dat');
    lenth = [3;7;15];
    par = [2;5;10];
end

switch q
    case 4
    v = xp(1:lenth(1), 1);
    u = xp(1:lenth(1), 2);
    a = par(1);
    case 8
    v = xp(1+lenth(1):lenth(1)+lenth(2),1);
    u = xp(1+lenth(1):lenth(1)+lenth(2),2);
    a = par(2);
    case 16
    v = xp(1+lenth(2)+lenth(1):sum(lenth),1);
    u = xp(1+lenth(2)+lenth(1):sum(lenth),2);
    a = par(3);
end