function y = normPlot(N,i,j)

figure()
plot(1:7,squeeze(N(i,j,:)), LineWidth=3)