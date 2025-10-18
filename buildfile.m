function plan = buildfile
import matlab.buildtool.tasks.CodeIssuesTask
import matlab.buildtool.tasks.TestTask
import matlab.buildtool.tasks.CleanTask
import matlab.buildtool.tasks.MexTask

% Create a plan from task functions
plan = buildplan(localfunctions);

% Add a task to identify code issues
plan("check") = CodeIssuesTask;

% Add a task to delete outputs and traces
plan("clean") = CleanTask;

% Add a task group to build MEX files
if ismac
    options = ["-ld_classic"];
else
    options = [];
end
plan("mex:DP") = MexTask("DP/DynamicProgrammingQ.c","mex");
plan("mex:DP2") = MexTask(["DP/DynamicProgrammingQ2.c", "DP/dp_grid.c", "DP/dp_nbhd.c"],"mex");
plan("mex:bcalcY") = MexTask("armadillo_cpp/bcalcY.cpp","mex",Options=options);
plan("mex:bcuL2norm2") = MexTask("armadillo_cpp/bcuL2norm2.cpp","mex",Options=options);
plan("mex:trapzCpp") = MexTask("armadillo_cpp/trapzCpp.cpp","mex",Options=options);
plan("mex:border_l2norm") = MexTask("armadillo_cpp/border_l2norm.cpp","mex",Options=options);
plan("mex:mcholC") = MexTask("minFunc/mex/mcholC.c", "mex");
plan("mex:lbfgsC") = MexTask("minFunc/mex/lbfgsC.c", "mex");
plan("mex:lbfgsAddC") = MexTask("minFunc/mex/lbfgsAddC.c", "mex");
plan("mex:lbfgsProdC") = MexTask("minFunc/mex/lbfgsProdC.c", "mex");
plan("mex:mlogit_warp") = MexTask(["mlogit_warp/mlogit_warp.c","mlogit_warp/mlogit_warp_grad.c", "mlogit_warp/misc_funcs.c"],"mex");
if ismac
    options = ["-ld_classic" "-llapack" "-lblas"];
elseif ispc
    libDir = fullfile(matlabroot,"/extern/lib/win64/microsoft");
    options = [sprintf('-L"%s"',libDir) "-ld_classic" "-lmwlapack" "-lmwblas"];
else
    options = ["-llapack" "-lblas"];
end
plan("mex:c_rlbfgs") = MexTask("armadillo_cpp/c_rlbfgs.cpp","mex",Options=options);


% Add a task to run tests
plan("test") = TestTask;

% Make the "archive" task the default task in the plan
plan.DefaultTasks = "package";

% Make the "archive" task dependent
plan("package").Dependencies = ["check" "clean" "mex" "test"];
end

function packageTask(~)
    % Define actions for packaging the toolbox
    disp("Packaging toolbox...");
    matlab.addons.toolbox.packageToolbox("fdasrvf.prj");
end
