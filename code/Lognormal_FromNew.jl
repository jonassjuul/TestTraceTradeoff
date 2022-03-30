
using DelimitedFiles

# Include my functions.....
include("Functions_Lognormal_FromNew.jl")


# Definitions
# -----------

#TimeMax = 40;     # Maximum days run.
InitialNumberOfInfected  = 100;      # Infectious at start
MaximumAllowedInfected = 100000; # How many people will we maximally get?

NumberOfExperiments  = 20; # Number of experiments

# Epidemiological details
AsymptomaticFractionOfInfected = 0.6;#0.3; # Fraction of infected that never get symptoms.

R0 = 2.0; #3//2.5 # Mean number of children in full period of infection.
OffspringDistribution = "poisson";
#OffspringDistribution = "geometric";

InfectiousProfile = "empirical";
#InfectiousProfile = "FlatSkewed";

MeanOfLognormal = getMeanOfLognormalDistribution();

# Societal details
WaitBeforeTestTaken  = 4;  # Number of days before test is taken
WaitBeforeTestResult  = 1; # Number of days before test result arrives after test is taken

# Test-and-trace details
ProbabilityChildIsTraced  = -0.02; #+34*0.02 // Fraction of children that are found through contact tracing
ProbabilityFalseNegativeTest = -0.02;

#--------------------
# Define directory where results will be saved
DirectoryToSaveResults = "Outputs/"

# Define Filename where results will be saved
FilenameToSaveResults = string("JULIA_TestSensitivity_Istart" ,InitialNumberOfInfected,"_Nexp",NumberOfExperiments,"_R0",R0,"_WaitBeforeTestTaken",WaitBeforeTestTaken,"_WaitBeforeTestResult",WaitBeforeTestResult, "_Asymptomatics",AsymptomaticFractionOfInfected,"_InfectiousProfile",InfectiousProfile,"_OffspringDistribution", OffspringDistribution,".txt");

# First list in filename where results will be saved specifies columns
FirstLineInFile = string("False negative test rate,","Tracing efficiency,","N_infected_done,","N_recovered,","Reff");
AppendLineToFile(string(DirectoryToSaveResults,FilenameToSaveResults),FirstLineInFile)
# Loop over different choices for 
#   1. Contact tracing efficiency (probability that a child is traced when parent gets tested positive.)
#   2. Test sensitivity

for TracingEfficiencyValueNumber = 1:50
    # Each time model is run for a new Tracing Efficiency Value, increase ProbabilityChildIsTraced
    global ProbabilityChildIsTraced += 0.02;
    global WaitBeforeTestTaken = WaitBeforeTestTaken;
    global WaitBeforeTestResult = WaitBeforeTestResult;

    # Each time Tracing Efficiency Value increases, reset ProbabilityFalseNegativeTest
    global ProbabilityFalseNegativeTest= -0.02;

    for TestSensitivityValueNumber =1:1#50
        # Each time model is run for a new Tracing Efficiency Value, increase ProbabilityFalseNegativeTest
        global ProbabilityFalseNegativeTest += 0.02;

        # Print progress.
        print("Currently simulating parameters:\t", "False neg:\t", ProbabilityFalseNegativeTest, "\tTrace efficiency:\t", ProbabilityChildIsTraced,"\n")

        


        # Do NumberOfExperiments runs for each parameter combination. Results will be average results over these experiments.
        
        # Define variables for averaged results
        RecoveredPeople_AveragedOverExperiments = 0;
        InfectedPeople_AveragedOverExperiments = 0;
        EffectiveReproductionNumber_AveragedOverExperiments = 0;

        # Do NumberOfExperiments runs for each parameter combination. 
        for ExperimentNumber =1:NumberOfExperiments
            #println("Experiment number\t",ExperimentNumber,"\tof:\t",NumberOfExperiments)
            # Define variables and vectors for each run.
            # -------

            # Integers
            NumberOfInfected = InitialNumberOfInfected +0; # Number of infected at the beginning of simulation.
            NumberOfRecovered = 0; # Number of people that recovered from disease.
            NumberOfPeopleDoneInfecting = 0; # Number of people that infected all that they will infect.
            # TO DO: CHeck difference between NumberOfRecovered and NumberOfPeopleDoneInfecting. Document this.


            # State arrays
            # TO DO: DEFINE THESE.
            StateOfNodes = zeros(MaximumAllowedInfected); # Array with 0, 1, 2, 3 on entry, corresponding to S, E, I, R.

            CountUpToStateChange = ones(MaximumAllowedInfected)*(-1) ;     # Array with Day on entry. Counts from 0. Node i changes state when Counter_goal[i] is reached.
            GoalOfCountDown = zeros(MaximumAllowedInfected) ;# Array with Day on entry. Negative if not infected,

            TestArrivalTimeOfNodes = ones(MaximumAllowedInfected)*(-9) ;     # Array with Day on entry. Negative if not waiting,
            ResultArrivalTimeOfNodes = ones(MaximumAllowedInfected)*(-9) ;   # Array with Day on entry. Negative if not waiting,
            TraceNodesChildren = zeros(MaximumAllowedInfected) ;             # Array with 0 or 1 on entry. 0 if not waiting to be traced,
            NodeCanTestPositive = zeros(MaximumAllowedInfected); # Array with 0 on entry if node cannot test positive. 1 if node can.

            Asymptomatic = floor.(Int,rand(MaximumAllowedInfected,1).+(AsymptomaticFractionOfInfected)); # Array with 0 or 1 on entry. 0 if normal, 1 if always asymptomatic.

            WhenInfectedWillInfectOthers = fill(Int[], MaximumAllowedInfected,1); # List at entry i contains days after infection when node i will infect other nodes.
            ListOfChildren = fill(Int[], MaximumAllowedInfected,1); # List at entry i contains nodes that node i infected. Used for contact tracing.

            # Infect a number of people at start of simulation
            StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers = getInitialConditionsOfSimulation(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,InitialNumberOfInfected,R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile);



            # Run model until noone is active anymore
            NodesStillActive=true;
            TimeStep = 0;
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
                StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren = TraceNode(StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren,Asymptomatic,MaximumAllowedInfected,WaitBeforeTestTaken,ProbabilityChildIsTraced);
            end
            # Having run simulation to end, count results in average over simulations.
            EffectiveReproductionNumber_AveragedOverExperiments += ((NumberOfInfected-InitialNumberOfInfected)/NumberOfRecovered)/NumberOfExperiments
            RecoveredPeople_AveragedOverExperiments += NumberOfRecovered/NumberOfExperiments
            InfectedPeople_AveragedOverExperiments += (NumberOfInfected-InitialNumberOfInfected)/NumberOfExperiments
            println("Experiment:\t",ExperimentNumber,"\tInfected:\t",NumberOfInfected,"\tRecovered:\t",NumberOfRecovered,"\tR0:\t",(NumberOfInfected-InitialNumberOfInfected)/NumberOfRecovered)
        end

        # Print averaged results to file.
        AveragedResultsToPrintToFile = string(ProbabilityFalseNegativeTest,",",ProbabilityChildIsTraced,",",InfectedPeople_AveragedOverExperiments,",",RecoveredPeople_AveragedOverExperiments,",",EffectiveReproductionNumber_AveragedOverExperiments);
        AppendLineToFile(string(DirectoryToSaveResults,FilenameToSaveResults),AveragedResultsToPrintToFile)


    end


end