suite = testsuite(fileparts(mfilename('fullpath')));
runner = testrunner;
results = run(runner,suite)
