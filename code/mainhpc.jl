
using DelimitedFiles
using Statistics
using Distributions
using Base.Threads
# Include my functions.....
include("Functions_main.jl")


# Definitions
# -----------

#TimeMax = 40;     # Maximum days run.
InitialNumberOfInfected  = 100;      # Infectious at start
MaximumAllowedInfected = 10000; # How many people will we maximally get?
N = MaximumAllowedInfected; # For theoretical calculations

NumberOfExperiments  = 500; # Number of experiments

# Epidemiological details
AsymptomaticFractionOfInfected = parse(Float64, ARGS[3]);# Fraction of infected that never get symptoms. 

R0 = parse(Float64, ARGS[4]); #3//2.5 # Mean number of children in full period of infection.

OffspringDistribution = ARGS[5];
# OffspringDistribution = "poisson";
# OffspringDistribution = "geometric";

InfectiousProfile = ARGS[6];
# InfectiousProfile = "empirical";
# InfectiousProfile = "FlatSkewed";
MeanOfLognormal = getMeanOfLognormalDistribution();

# Societal details
WaitBeforeTestTaken  = parse(Float64, ARGS[1]);  # Number of days before test is taken
WaitBeforeTestResult  = parse(Float64, ARGS[2]); # Number of days before test result arrives after test is taken

# Test-and-trace details
ProbabilityChildIsTraced  = -0.02; #+34*0.02 // Fraction of children that are found through contact tracing
ProbabilityFalseNegativeTest = -0.02;
linspace = 51;
#--------------------
# Define directory where results will be saved
DirectoryToSaveResults = "code/OutputsHPC/";

# Define Filename where results will be saved
FilenameToSaveResults = string("JULIA_TestSensitivity_Istart" ,InitialNumberOfInfected,"_Nexp",NumberOfExperiments,"_R0",R0,"_WaitBeforeTestTaken",Int(WaitBeforeTestTaken),"_WaitBeforeTestResult",Int(WaitBeforeTestResult), "_Asymptomatics",AsymptomaticFractionOfInfected,"_InfectiousProfile",InfectiousProfile,"_OffspringDistribution", OffspringDistribution,".txt");

# First list in filename where results will be saved specifies columns
# FirstLineInFile = string("False negative test rate,","Tracing efficiency,","N_infected_done,","N_recovered,","ReffMean,","ReffStd,","ReffTheoretical,","N_traced");

#Check if file already exists and handle accordingly
# if isfile(string(DirectoryToSaveResults,FilenameToSaveResults))
#     println("ERROR: File already exists: ", string(DirectoryToSaveResults,FilenameToSaveResults))
#     println("Please remove the existing file or change the filename to avoid overwriting data.")
#     error("Execution stopped to prevent overwriting existing file.")
# end

AppendLineToFile(string(DirectoryToSaveResults,FilenameToSaveResults),FirstLineInFile);

# Loop over different choices for 
#   1. Contact tracing efficiency (probability that a child is traced when parent gets tested positive.)
#   2. Test sensitivity

if WaitBeforeTestTaken + WaitBeforeTestResult == 0 #only run with one test sensitivity value if there is a delay (slow test)
    linspace = 26
    ProbabilityChildIsTraced = 0.5000000000000001
end

elapsed_time = @elapsed for TracingEfficiencyValueNumber = 1:linspace
    # Each time model is run for a new Tracing Efficiency Value, increase ProbabilityChildIsTraced
    global ProbabilityChildIsTraced += 0.02;
    global WaitBeforeTestTaken = WaitBeforeTestTaken;
    global WaitBeforeTestResult = WaitBeforeTestResult;

    # Each time Tracing Efficiency Value increases, reset ProbabilityFalseNegativeTest
    global ProbabilityFalseNegativeTest = -0.02;

    if WaitBeforeTestTaken + WaitBeforeTestResult > 0 #only run with one test sensitivity value if there is a delay (slow test)
        testsentivitylinspace = 1
    else
        testsentivitylinspace = linspace
    end

    for TestSensitivityValueNumber = 1:51
        # Each time model is run for a new Tracing Efficiency Value, increase ProbabilityFalseNegativeTest
        global ProbabilityFalseNegativeTest += 0.02;
        # Print progress.
        print("\nCurrently simulating parameters:\t", "False neg:\t", ProbabilityFalseNegativeTest, "\tTrace efficiency:\t", ProbabilityChildIsTraced,"\n")

        # Do NumberOfExperiments runs for each parameter combination. Results will be average results over these experiments.
        
        # Define variables for averaged results
        Recovered_each = zeros(Float64, NumberOfExperiments)
        Infected_each = zeros(Float64, NumberOfExperiments)
        Reff_each = zeros(Float64, NumberOfExperiments)
        PropTraced_each = zeros(Float64, NumberOfExperiments)
        GoalOfCountDown_traced_each = [Any[] for _ in 1:NumberOfExperiments]
        timetraced_traced_each = [Any[] for _ in 1:NumberOfExperiments]
        GoalOfCountDown_untraced_each = [Any[] for _ in 1:NumberOfExperiments] #store conditioned infection period lengths for untraced nodes for estimating the theoretical Reff
        # Do NumberOfExperiments runs for each parameter combination. 
        for ExperimentNumber = 1:NumberOfExperiments
            #println("Experiment number\t",ExperimentNumber,"\tof:\t",NumberOfExperiments)
            # Define variables and vectors for each run.
            # -------

            # Integers
            NumberOfInfected = InitialNumberOfInfected +0; # Number of infected at the beginning of simulation.
            NumberOfRecovered = 0; # Number of people that recovered from disease.
            NumberOfPeopleDoneInfecting = 0; # Number of people that infected all that they will infect.
            sumInfectiontimeGivenTracing = 0.0
            # TO DO: CHeck difference between NumberOfRecovered and NumberOfPeopleDoneInfecting. Document this.

            # State arrays
            # TO DO: DEFINE THESE.
            StateOfNodes = zeros(MaximumAllowedInfected); # Array with 0, 1, 2, 3 on entry, corresponding to S, E, I, R.

            CountUpToStateChange = ones(MaximumAllowedInfected)*(-1) ;     # Array with Day on entry. Counts from 0. Node i changes state when Counter_goal[i] is reached.
            GoalOfCountDown = zeros(MaximumAllowedInfected) ;# Array with Day on entry. Negative if not infected,

            TestArrivalTimeOfNodes = ones(MaximumAllowedInfected)*(-9) ;     # Array with Day on entry. Negative if not waiting,
            ResultArrivalTimeOfNodes = ones(MaximumAllowedInfected)*(-9) ;   # Array with Day on entry. Negative if not waiting,
            TraceNodesChildren = zeros(MaximumAllowedInfected);         # Array with 0 or 1 on entry. 0 if not waiting to be traced,
            NodeCanTestPositive = zeros(MaximumAllowedInfected); # Array with 0 on entry if node cannot test positive. 1 if node can.
            tracedNodes = zeros(MaximumAllowedInfected) #0 if node is never traced, 1 if node is traced, 2 if node is traced to late
            Asymptomatic = floor.(Int,rand(MaximumAllowedInfected,1).+(AsymptomaticFractionOfInfected)); # Array with 0 or 1 on entry. 0 if normal, 1 if always asymptomatic.

            WhenInfectedWillInfectOthers = fill(Int[], MaximumAllowedInfected,1); # List at entry i contains days after infection when node i will infect other nodes.
            ListOfChildren = fill(Int[], MaximumAllowedInfected,1); # List at entry i contains nodes that node i infected. Used for contact tracing.

            # Add this after the simulation completes to count entries
            

            # Infect a number of people at start of simulation
            StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers = getInitialConditionsOfSimulation(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,InitialNumberOfInfected,R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile);



            # Run model until none is active anymore
            NodesStillActive=true;
            TimeStep = 0;
            WaitBeforeTestResult_local = WaitBeforeTestResult
            GoalOfCountDown_traced_local = Any[]
            timetraced_traced_local = Any[]
            GoalOfCountDown_untraced_local = Any[]
            while NodesStillActive==true
                # Advance Time 1 step
                TimeStep +=1;
                # Advance all infected and all waiting 1 time step.
                StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,NodeCanTestPositive,TraceNodesChildren,WaitBeforeTestResult,NumberOfRecovered,FoundNoInfectiousOrExposedNode = AdvanceInfectedOneTimestep(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,NodeCanTestPositive,TraceNodesChildren,MaximumAllowedInfected,R0,MeanOfLognormal,WaitBeforeTestResult,NumberOfRecovered,ProbabilityFalseNegativeTest,OffspringDistribution,InfectiousProfile)

                # If no nodes are infectiuos or exposed, stop simulation.
                if FoundNoInfectiousOrExposedNode == true
                    NodesStillActive = false;
                    break
                end

                # Infect all children that are due to get infected this time step.
                StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,ListOfChildren,NumberOfInfected=InfectNodesOnThisTimestep(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,MaximumAllowedInfected,NumberOfInfected,R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile);

                # Trace nodes that should get traced this time step and test nodes that get symptoms.
                StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren, tracedNodes, GoalOfCountDown_traced_local, timetraced_traced_local, GoalOfCountDown_untraced_local = TraceNode(StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren,Asymptomatic,MaximumAllowedInfected,WaitBeforeTestTaken,ProbabilityChildIsTraced,tracedNodes,GoalOfCountDown_traced_local,timetraced_traced_local, GoalOfCountDown_untraced_local);
            end

             #count proportion of traced nodes
            if NumberOfInfected < MaximumAllowedInfected
                PropTraced_each[ExperimentNumber] = sum(tracedNodes .== 1)/(NumberOfInfected-InitialNumberOfInfected)
            else
                PropTraced_each[ExperimentNumber] = sum(tracedNodes .== 1)/(MaximumAllowedInfected-InitialNumberOfInfected)
            end

            # Having run simulation to end, count results in average over simulations.
            Reff_each[ExperimentNumber] = (NumberOfInfected-InitialNumberOfInfected)/NumberOfRecovered
            Recovered_each[ExperimentNumber] = NumberOfRecovered
            Infected_each[ExperimentNumber] = (NumberOfInfected-InitialNumberOfInfected)
            GoalOfCountDown_traced_each[ExperimentNumber] = GoalOfCountDown_traced_local
            timetraced_traced_each[ExperimentNumber] = timetraced_traced_local
            GoalOfCountDown_untraced_each[ExperimentNumber] = GoalOfCountDown_untraced_local
        end

        RecoveredPeople_AveragedOverExperiments = mean(Recovered_each)
        InfectedPeople_AveragedOverExperiments = mean(Infected_each)
        EffectiveReproduction = Reff_each
        avgProportionTraced = mean(PropTraced_each)
        GoalOfCountDown_traced = vcat(GoalOfCountDown_traced_each...)
        timetraced_traced = vcat(timetraced_traced_each...)
        GoalOfCountDown_untraced = vcat(GoalOfCountDown_untraced_each...)


        τ = Int(WaitBeforeTestTaken + WaitBeforeTestResult)
        I1untrace = 0.0
        I2untrace = 0.0

        count = 0
        Nuntrace = length(GoalOfCountDown_untraced)
        randomVariable = 0
        for test in 1:Nuntrace
            t_half = GoalOfCountDown_untraced[test] /2 + 1
            InfectiousnessDistribution = getInfectiousnessDistribution(GoalOfCountDown_untraced[test], InfectiousProfile);
            R0OfNode = R0 * GoalOfCountDown_untraced[test] / (2*MeanOfLognormal);
            NumberOfChildrenToDraw = drawNumberOfChildren(R0OfNode,OffspringDistribution)
            L = length(InfectiousnessDistribution);

            # Convert to integer indices
            b1 = floor(Int, t_half);
            b2 = b1 + τ;

            if b2 <= L
                I1untrace += sum(InfectiousnessDistribution[b2+1:end])*NumberOfChildrenToDraw;
                I2untrace += sum(InfectiousnessDistribution[b1+1:b2])*NumberOfChildrenToDraw;
            else
                I1untrace += 0;
                I2untrace += sum(InfectiousnessDistribution[b1+1:end])*NumberOfChildrenToDraw;
                count += 1;
            end
        end

        theoreticalMean_NoTracingTerm = (1-AsymptomaticFractionOfInfected)*((1-ProbabilityFalseNegativeTest)*I1untrace/Nuntrace+I2untrace/Nuntrace)
        
        I1trace = 0.0
        I2trace = 0.0
        Ntrace = length(timetraced_traced)
        if Ntrace != 0
            for i in 1:Ntrace
                InfectiousnessDistribution = getInfectiousnessDistribution(GoalOfCountDown_traced[i], InfectiousProfile);
                R0OfNode = R0 * GoalOfCountDown_traced[i] / (2*MeanOfLognormal);
                NumberOfChildrenToDraw = drawNumberOfChildren(R0OfNode,OffspringDistribution)
                L = length(InfectiousnessDistribution);
                b1 = Int(timetraced_traced[i]);
                b2 = b1 + τ;

                if b2 <= L
                    I1trace += sum(InfectiousnessDistribution[b2+1:end])*NumberOfChildrenToDraw;
                    I2trace += sum(InfectiousnessDistribution[b1+1:b2])*NumberOfChildrenToDraw;
                else
                    I1trace += 0;
                    I2trace += sum(InfectiousnessDistribution[b1+1:end])*NumberOfChildrenToDraw;
                end
            end

            theoreticalMean_TracingTerm = ((1-ProbabilityFalseNegativeTest)*I1trace/Ntrace+I2trace/Ntrace)
        else
            theoreticalMean_TracingTerm = 0.0
        end
        
        if Ntrace != 0
            pTraceEff = avgProportionTraced 
        else
            pTraceEff = 0.0
        end

        theoreticalMean = R0 - (1-pTraceEff)*theoreticalMean_NoTracingTerm - pTraceEff*theoreticalMean_TracingTerm
        print("Theoretical Reff:\t", theoreticalMean, "\n")
        print("Empirical Reff:\t", mean(EffectiveReproduction), "\n")
        # Print averaged results to file.
        AveragedResultsToPrintToFile = string(ProbabilityFalseNegativeTest,",",ProbabilityChildIsTraced,",",InfectedPeople_AveragedOverExperiments,",",RecoveredPeople_AveragedOverExperiments,",",mean(EffectiveReproduction), ",", std(EffectiveReproduction)/sqrt(NumberOfExperiments), ",", theoreticalMean, ",", Ntrace/NumberOfExperiments);
        AppendLineToFile(string(DirectoryToSaveResults,FilenameToSaveResults),AveragedResultsToPrintToFile)
  
    end
end
print("\n Total elapsed time:\t", elapsed_time, " seconds.\n")


