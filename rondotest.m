% create parallel pool
%
if (strfind(version('-date'), '2013')) % rondo lives in 2013
    matlabpool;
    fprintf('num of 2013 parallel workers = %d\n', matlabpool('size'));
else
    poolobj = gcp;
    fprintf('num of parallel workers = %d\n', poolobj.NumWorkers);
end


a = 1:100;

parfor i=a
    i
end


init_par = [0.35    0.3, ...     % focal, low emph 
            0.6     0.4, ...     % focal, high emph 
            0.8     0.75, ...    % nonfocal, low emph 
            0.9     0.83, ...    % nonfocal, high emph
            0.4]; % gamma * 10^-3

error = fitme(init_par);

error

fprintf('error = %f\n', error);
