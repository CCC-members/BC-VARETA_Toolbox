% define MATLAB matrices
pyenv('Version','/usr/local/bin/python3.9',ExecutionMode="OutOfProcess");

A1=[1 2;
    3 4];
A2=[2 3;
    1 0];
 
% check the conversion types
% https://www.mathworks.com/help/matlab/matlab_external/passing-data-to-python.html
% convert matrices to Python objects
A1converted = py.numpy.array(A1);
A2converted = py.numpy.array(A2);
% parameter to be sent
parameter=py.int(2);