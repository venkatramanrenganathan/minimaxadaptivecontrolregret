deltas = 0.5:0.02:0.99;
risks = norminv(1-deltas);
drrisks = zeros(size(deltas,2),1);
for i = 1:size(deltas,2)
    drrisks(i, 1) = sqrt(deltas(1,i)/(1 - deltas(1,i))); 
end
plot(deltas, risks, '-b');
hold on;
plot(deltas, drrisks, '-r');