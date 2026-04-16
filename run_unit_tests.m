clear; clc;
% runMyTests.m
import matlab.unittest.TestSuite;

% 1. Add current folder and Unit_tests to path
addpath(pwd);
addpath(fullfile(pwd, 'Unit_Test'));

% % 2. Create the suite from the class name
% suite = TestSuite.fromClass(?ECEF_ECI_Test);
% 
% % 3. Run and view results
% result = run(suite);
% display(result);

suite = TestSuite.fromClass(?Sun_Moon_Vector_Tests);
result = run(suite);
display(result);
