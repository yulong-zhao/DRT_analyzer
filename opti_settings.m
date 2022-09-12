opti.fmincon = optimoptions('fmincon', 'MaxFunctionEvaluations', 300e3);
opti.fmincon.Display = 'final';
opti.fmincon.MaxIterations = 300e3;
opti.fmincon.UseParallel = true;

opti.fminsearch = optimset();
opti.fminsearch.Display = 'final'; %iter-detailed
opti.fminsearch.MaxFunEvals = 200e3;
opti.fminsearch.MaxIter = 500e3;
opti.fmincon.UseParallel = true;


opti.patternsearch = optimoptions('patternsearch');
opti.patternsearch.Display = 'iter';
opti.patternsearch.MaxIterations = 300e3;
opti.patternsearch.MaxFunctionEvaluations = 300e3;
opti.patternsearch.UseParallel = true;


opti.particleswarm = optimoptions('particleswarm', 'SwarmSize', round(10e3));
opti.particleswarm.HybridFcn = {@fmincon,opti.fmincon};
opti.particleswarm.Display = 'iter';
opti.particleswarm.InitialSwarmMatrix = PSO_InitialSwarmMatrix;
opti.particleswarm.UseParallel = true;
opti.particleswarm.MaxTime = appStruct.other.MaxTime;

opti.simulannealbnd = optimoptions('simulannealbnd');
opti.simulannealbnd.Display = 'iter';


opti.lsqnonlin = optimoptions('lsqnonlin');
opti.lsqnonlin.Display = 'final';
opti.lsqnonlin.MaxFunctionEvaluations = 300e3;
opti.lsqnonlin.MaxIterations = 100e3;
opti.lsqnonlin.UseParallel = true;

%opti.lsqnonlin.Algorithm = 'levenberg-marquardt';


% --- custom p





opti.ga = optimoptions('ga','InitialPopulationMatrix', PSO_InitialSwarmMatrix, 'PopulationSize', 3000); 
opti.ga.Display = 'iter';

opti.ga.MaxTime = appStruct.other.MaxTime;
opti.ga.UseParallel = true;
opti.ga.HybridFcn = {@fmincon,opti.fmincon};
