clear; clc;

import matlab.unittest.TestSuite;

% 1. Setup Paths
addpath(genpath(pwd)); % Adds everything including Unit_Test folder

% 2. Create a Suite from the entire Unit_Test folder
% This automatically finds all files ending in "Test.m" or containing test classes
suite = TestSuite.fromFolder(fullfile(pwd, 'Unit_Test'));

% 3. Create a TestRunner for better reporting
runner = matlab.unittest.TestRunner.withTextOutput;

% Optional: Add a PDF or HTML report plugin
% reportFile = 'TestResults.pdf';
% runner.addPlugin(matlab.unittest.plugins.TestReportPlugin.producingPDF(reportFile));

% 4. Run the suite
results = runner.run(suite);

% 5. Format the output for a quick "Pass/Fail" overview
table(results)