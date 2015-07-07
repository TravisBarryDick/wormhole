centers = 3.5 * [-1 -1; 0 -1; 1 -1 ;  -1 0; 0 0; 1 0; -1 1; 0 1; 1 1;];

dataset = [];
f  = 0.25;
n = 80;
for (i=centers')
    %display(i')
    dataset = [dataset ; mvnrnd(i', f * eye(2), n)];
end
dataset  = dataset(randperm(size(dataset, 1)),:);



[y, X] = libsvmread('/Users/krishnapillutla/Dropbox/Work/Acad/DistributedLearning/code/c++/sample_out');
X1 = full(X);
pts = cell(9, 1);
for (i=1:9)
   pts{i} = X1((y == i-1), :); 
   size(pts{i})
end
