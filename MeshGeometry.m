
A = importdata('cube.txt');

scatter3(A(:,1),A(:,2),A(:,3),10,A(:,3))
axis([-20 20 -20 20 -20 20])
colormap(jet);
colorbar;














