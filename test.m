hold on
for i=1:length(lambda)-1
    rms(NA(v(:,i),k,a,b,c)-findRHS(x,lambda(i)))
    plot(x,NA(v(:,i),k,a,b,c)-findRHS(x,lambda(i)))
end
hold off


