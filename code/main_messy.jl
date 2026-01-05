
using DelimitedFiles
using Statistics
# Include my functions.....
include("Functions_main_messy.jl")


# Definitions
# -----------

#TimeMax = 40;     # Maximum days run.
InitialNumberOfInfected  = 100;      # Infectious at start
MaximumAllowedInfected = 100000; # How many people will we maximally get?
N = MaximumAllowedInfected; # For theoretical calculations

NumberOfExperiments  = 50; # Number of experiments

# Epidemiological details
AsymptomaticFractionOfInfected = 0.3;# Fraction of infected that never get symptoms. 

R0 = 2.0; #3//2.5 # Mean number of children in full period of infection.
OffspringDistribution = "poisson";
# OffspringDistribution = "geometric";

InfectiousProfile = "empirical";
#InfectiousProfile = "FlatSkewed";
MeanOfLognormal = getMeanOfLognormalDistribution();

# Societal details
WaitBeforeTestTaken  = 0;  # Number of days before test is taken
WaitBeforeTestResult  = 0; # Number of days before test result arrives after test is taken

# Test-and-trace details
ProbabilityChildIsTraced  = 0.78; #+34*0.02 // Fraction of children that are found through contact tracing
ProbabilityFalseNegativeTest = -0.02;
linspace = 51;
empiricalMean = zeros(linspace, linspace);
theoreticalMeanArray = zeros(linspace, linspace);
averageDayofInfectionMatrix = zeros(linspace, linspace);
averageDayofInfectionTracingMatrix = zeros(linspace, linspace);
averagenumbertracedMatrix = zeros(linspace,linspace);
avgProportionTracedMatrix = zeros(linspace,linspace);
avgPandemicEnded = zeros(linspace,linspace);
PTraceValidMatrix = zeros(linspace,linspace);
avgProportionTracedTheoretical = zeros(linspace,linspace);
#--------------------
# Define directory where results will be saved
DirectoryToSaveResults = "code/OutputspfalseDistribution/";

# Define Filename where results will be saved
FilenameToSaveResults = string("JULIA_TestSensitivity_Istart" ,InitialNumberOfInfected,"_Nexp",NumberOfExperiments,"_R0",R0,"_WaitBeforeTestTaken",WaitBeforeTestTaken,"_WaitBeforeTestResult",WaitBeforeTestResult, "_Asymptomatics",AsymptomaticFractionOfInfected,"_InfectiousProfile",InfectiousProfile,"_OffspringDistribution", OffspringDistribution,".txt");

# First list in filename where results will be saved specifies columns
FirstLineInFile = string("False negative test rate,","Tracing efficiency,","N_infected_done,","N_recovered,","ReffMean,","ReffStd,","ReffTheoretical,","N_traced");
AppendLineToFile(string(DirectoryToSaveResults,FilenameToSaveResults),FirstLineInFile)

# Loop over different choices for 
#   1. Contact tracing efficiency (probability that a child is traced when parent gets tested positive.)
#   2. Test sensitivity

elapsed_time = @elapsed for TracingEfficiencyValueNumber = 1:21
    # Each time model is run for a new Tracing Efficiency Value, increase ProbabilityChildIsTraced
    global ProbabilityChildIsTraced += 0.02;
    global WaitBeforeTestTaken = WaitBeforeTestTaken;
    global WaitBeforeTestResult = WaitBeforeTestResult;

    # Each time Tracing Efficiency Value increases, reset ProbabilityFalseNegativeTest
    global ProbabilityFalseNegativeTest = -0.02;

    for TestSensitivityValueNumber =1:linspace
        # Each time model is run for a new Tracing Efficiency Value, increase ProbabilityFalseNegativeTest
        global ProbabilityFalseNegativeTest += 0.02;
        # Print progress.
        print("\nCurrently simulating parameters:\t", "False neg:\t", ProbabilityFalseNegativeTest, "\tTrace efficiency:\t", ProbabilityChildIsTraced,"\n")

        # Do NumberOfExperiments runs for each parameter combination. Results will be average results over these experiments.
        
        # Define variables for averaged results
        RecoveredPeople_AveragedOverExperiments = 0;
        InfectedPeople_AveragedOverExperiments = 0;
        EffectiveReproduction = [];
        averageDayofInfection = 0;
        averageDayofInfectionGivenTracing = 0;
        endedPandemics = 0; 
        avgnumbertraced = 0;
        avgProportionTraced = 0;
        GoalOfCountDown_traced = [];
        timetraced_traced = [];
        GoalOfCountDown_untraced = [];
        nonvalidtracings = 0;
        PTraceValidArrayexp = 0.0;
        nExpWithTracings = 0;
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
            pTestDistribution = fill(Float64[], MaximumAllowedInfected,1);
            # Add this after the simulation completes to count entries
            

            # Infect a number of people at start of simulation
            StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers, pTestDistribution = getInitialConditionsOfSimulation(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,InitialNumberOfInfected,R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile, pTestDistribution, ProbabilityFalseNegativeTest);



            # Run model until none is active anymore
            NodesStillActive=true;
            global TimeStep = 0;
            while NodesStillActive==true
                # Advance Time 1 step
                TimeStep +=1;
                # Advance all infected and all waiting 1 time step.
                StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,NodeCanTestPositive,TraceNodesChildren,WaitBeforeTestResult,NumberOfRecovered,FoundNoInfectiousOrExposedNode = AdvanceInfectedOneTimestep(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,NodeCanTestPositive,TraceNodesChildren,MaximumAllowedInfected,R0,MeanOfLognormal,WaitBeforeTestResult,NumberOfRecovered,ProbabilityFalseNegativeTest,OffspringDistribution,InfectiousProfile, pTestDistribution)

                # If no nodes are infectiuos or exposed, stop simulation.
                if FoundNoInfectiousOrExposedNode == true
                    NodesStillActive = false;
                    break
                end

                # Infect all children that are due to get infected this time step.
                StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,ListOfChildren,NumberOfInfected, pTestDistribution=InfectNodesOnThisTimestep(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,MaximumAllowedInfected,NumberOfInfected,R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile, pTestDistribution, ProbabilityFalseNegativeTest);

                # Trace nodes that should get traced this time step and test nodes that get symptoms.
                StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren, sumInfectiontimeGivenTracing, tracedNodes, GoalOfCountDown_traced, timetraced_traced, nonvalidtracings, GoalOfCountDown_untraced = TraceNode(StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,WhenInfectedWillInfectOthers,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren,Asymptomatic,MaximumAllowedInfected,WaitBeforeTestTaken,ProbabilityChildIsTraced,sumInfectiontimeGivenTracing,tracedNodes,GoalOfCountDown_traced,timetraced_traced, nonvalidtracings, GoalOfCountDown_untraced);
            end

            if NumberOfInfected < MaximumAllowedInfected
                endedPandemics += 1
                avgProportionTraced += sum(tracedNodes .== 1)/(NumberOfInfected-InitialNumberOfInfected)/NumberOfExperiments   
            else
                avgProportionTraced += sum(tracedNodes .== 1)/(MaximumAllowedInfected-InitialNumberOfInfected)/NumberOfExperiments
            end

            if sum(tracedNodes .== 1) + nonvalidtracings != 0
                PTraceValidArrayexp += sum(tracedNodes .== 1)/(sum(tracedNodes .== 1)+nonvalidtracings)
                nExpWithTracings += 1
            end

            # Having run simulation to end, count results in average over simulations.
            push!(EffectiveReproduction, (NumberOfInfected-InitialNumberOfInfected)/NumberOfRecovered)
            RecoveredPeople_AveragedOverExperiments += NumberOfRecovered/NumberOfExperiments
            InfectedPeople_AveragedOverExperiments += (NumberOfInfected-InitialNumberOfInfected)/NumberOfExperiments
            averageDayofInfection += sum(sum.(WhenInfectedWillInfectOthers))/sum(length.(WhenInfectedWillInfectOthers))/NumberOfExperiments
            averageDayofInfectionGivenTracing += sumInfectiontimeGivenTracing/sum(tracedNodes .== 1)/NumberOfExperiments
            avgnumbertraced += sum(tracedNodes .== 1)/NumberOfExperiments
            #avgProportionTraced += totalChildrenTraced/(NumberOfInfected-InitialNumberOfInfected)/NumberOfExperiments
            #print("total Children traced:\t", totalChildrenTraced,"\tsum Infection time given tracing:\t", sumInfectiontimeGivenTracing,"\n")
            #println("Experiment:\t",ExperimentNumber,"\tInfected:\t",NumberOfInfected,"\tRecovered:\t",NumberOfRecovered,"\tR0:\t",(NumberOfInfected-InitialNumberOfInfected)/NumberOfRecovered)
            #print("Experiment ", ExperimentNumber
            #print(NumberOfRecovered," number of recovered\n", NumberOfInfected," number of infected\n")

        end
        #print("\n Ended pandemics:\t",endedPandemics," out of ",NumberOfExperiments,"\n")
        #print("Average day of infection:\t",averageDayofInfection,"\n")
        #print("Average day of infection given tracing:\t",averageDayofInfectionGivenTracing,"\n")
        # averageDayofInfectionMatrix[TracingEfficiencyValueNumber, TestSensitivityValueNumber] = averageDayofInfection
        # averageDayofInfectionTracingMatrix[TracingEfficiencyValueNumber, TestSensitivityValueNumber] = averageDayofInfectionGivenTracing
        # averagenumbertracedMatrix[TracingEfficiencyValueNumber, TestSensitivityValueNumber] = avgnumbertraced
        # avgProportionTracedMatrix[TracingEfficiencyValueNumber, TestSensitivityValueNumber] = avgProportionTraced
        # avgPandemicEnded[TracingEfficiencyValueNumber, TestSensitivityValueNumber] = endedPandemics/NumberOfExperiments
        # PTraceValidMatrix[TracingEfficiencyValueNumber, TestSensitivityValueNumber] = PTraceValidArrayexp/nExpWithTracings
        
        τ = WaitBeforeTestTaken + WaitBeforeTestResult 
        I1 = 0.0
        I2 = 0.0

        count = 0
        Nuntrace = length(GoalOfCountDown_untraced)
        randomVariable = 0
        for test in 1:Nuntrace
            t_half = GoalOfCountDown_untraced[test] /2 + 1
            InfectiousnessDistribution = getInfectiousnessDistribution(GoalOfCountDown_untraced[test], InfectiousProfile);
            R0OfNode = R0 * GoalOfCountDown_untraced[test] / (2*MeanOfLognormal);
            NumberOfChildrenToDraw = drawNumberOfChildren(R0OfNode,"poisson")
            L = length(InfectiousnessDistribution);

            # Convert to integer indices
            b1 = floor(Int, t_half);
            b2 = b1 + τ;

            if b2 <= L
                I1 += sum(InfectiousnessDistribution[b2+1:end])*NumberOfChildrenToDraw;
                I2 += sum(InfectiousnessDistribution[b1+1:b2])*NumberOfChildrenToDraw;
            else
                I1 += 0;
                I2 += sum(InfectiousnessDistribution[b1+1:end])*NumberOfChildrenToDraw;
                count += 1;
            end
        end
        print("Average untraced variable:\t", sum(GoalOfCountDown_untraced)/Nuntrace, "\n")



        theoreticalMean_NoTracingTerm = (1-AsymptomaticFractionOfInfected)*((1-ProbabilityFalseNegativeTest)*I1/Nuntrace+I2/Nuntrace)
        
        I1trace = 0.0
        I2trace = 0.0
        Ntrace = length(timetraced_traced)
        if Ntrace != 0
            for i in 1:Ntrace
                InfectiousnessDistribution = getInfectiousnessDistribution(GoalOfCountDown_traced[i], InfectiousProfile);
                R0OfNode = R0 * GoalOfCountDown_traced[i] / (2*MeanOfLognormal);
                NumberOfChildrenToDraw = drawNumberOfChildren(R0OfNode,"poisson")
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

        pTraceEff = (1-AsymptomaticFractionOfInfected)*ProbabilityChildIsTraced*(1-ProbabilityFalseNegativeTest)/(1-AsymptomaticFractionOfInfected*ProbabilityChildIsTraced*(1-ProbabilityFalseNegativeTest))
        avgProportionTracedTheoretical[TracingEfficiencyValueNumber, TestSensitivityValueNumber] = PTraceValidArrayexp*(1-AsymptomaticFractionOfInfected)*ProbabilityChildIsTraced*(1-ProbabilityFalseNegativeTest)/(1-AsymptomaticFractionOfInfected*ProbabilityChildIsTraced*(1-ProbabilityFalseNegativeTest)*PTraceValidArrayexp)
        if Ntrace != 0
            # print("valid traces proportion version 2:\t", PTraceValidArrayexp, "\n")
            pTraceEff = (PTraceValidArrayexp/nExpWithTracings)*(1-AsymptomaticFractionOfInfected)*ProbabilityChildIsTraced*(1-ProbabilityFalseNegativeTest)/(1-AsymptomaticFractionOfInfected*ProbabilityChildIsTraced*(1-ProbabilityFalseNegativeTest)*(PTraceValidArrayexp/nExpWithTracings))
            #empirical value
            # print("Adjusted pTraceEff:\t", pTraceEff, "\t Empiric value in comparison: ",avgProportionTraced,  "\n")
            # print("Difference theoretical-empirical proportion traced: \t", avgProportionTraced - pTraceEff, "\n")
            pTraceEff = avgProportionTraced 
        end

        theoreticalMean = (R0 - (1-pTraceEff)*theoreticalMean_NoTracingTerm - pTraceEff*theoreticalMean_TracingTerm)

        print("Theoretical Mean: ", theoreticalMean, "\n")
        print("Empirical Mean: ", mean(EffectiveReproduction),"+-", std(EffectiveReproduction),"\n")
        
        #TheoreticalMeanv2 = R0*AsymptomaticFractionOfInfected + R0*(1-AsymptomaticFractionOfInfected)*((1-ProbabilityFalseNegativeTest)*I3/N+ProbabilityFalseNegativeTest*(1-I2/N))
        # Print averaged results to file.
        AveragedResultsToPrintToFile = string(ProbabilityFalseNegativeTest,",",ProbabilityChildIsTraced,",",InfectedPeople_AveragedOverExperiments,",",RecoveredPeople_AveragedOverExperiments,",",mean(EffectiveReproduction), ",", std(EffectiveReproduction), ",", theoreticalMean, ",", Ntrace);
        AppendLineToFile(string(DirectoryToSaveResults,FilenameToSaveResults),AveragedResultsToPrintToFile)
  
        # print("\nTheoretical Reff:\t",theoreticalMean,"\n");

        # print("Count:\t",count,"\n Integral 1: ", I1/N, "\t Integral 2: ", I2/N, "\n");
        # print("\nTheoretical Reff v2:\t",TheoreticalMeanv2,"\n");
        # print("\nEffective Reff after simulation:\t",mean(EffectiveReproduction),"\n");
        #theoreticalMeanv2[TestSensitivityValueNumber] = TheoreticalMeanv2


    end
end
print("\n Total elapsed time:\t", elapsed_time, " seconds.\n")

writedlm(string(DirectoryToSaveResults, "empiricalMean_wait$(WaitBeforeTestResult + WaitBeforeTestTaken)NumExp$(NumberOfExperiments)withSampledPefftraced.txt"), empiricalMean)
writedlm(string(DirectoryToSaveResults, "theoreticalMean_wait$(WaitBeforeTestResult + WaitBeforeTestTaken)NumExp$(NumberOfExperiments)withSampledPefftraced.txt"), theoreticalMeanArray)

# # Save matrices to files
# writedlm(string(DirectoryToSaveResults, "averageDayofInfectionMatrix_wait$(WaitBeforeTestResult + WaitBeforeTestTaken).txt"), averageDayofInfectionMatrix)
# writedlm(string(DirectoryToSaveResults, "averageDayofInfectionTracingMatrix_wait$(WaitBeforeTestResult + WaitBeforeTestTaken).txt"), averageDayofInfectionTracingMatrix)
# writedlm(string(DirectoryToSaveResults,"averagenumbertracedMatrix_wait$(WaitBeforeTestResult + WaitBeforeTestTaken).txt"),averagenumbertracedMatrix)
# writedlm(string(DirectoryToSaveResults,"averageproportiontracedMatrix_wait$(WaitBeforeTestResult + WaitBeforeTestTaken).txt"),avgProportionTracedMatrix)
# writedlm(string(DirectoryToSaveResults,"averagepandemicended_wait$(WaitBeforeTestResult + WaitBeforeTestTaken).txt"),avgPandemicEnded)



using Plots

# Create x and y axis labels
tracing_efficiency_labels = 0:0.1:1.0
test_sensitivity_labels = 0:0.1:1.0

# Plot heatmap for average day of infection matrix
p1 = heatmap(test_sensitivity_labels, tracing_efficiency_labels, averageDayofInfectionMatrix,
    xlabel="Test Sensitivity (1 - False Negative Rate)",
    ylabel="Tracing Efficiency", 
    title="Average Day of Infection without Tracing \n",
    color=:viridis)

# Plot heatmap for average day of infection with tracing matrix
p2 = heatmap(test_sensitivity_labels, tracing_efficiency_labels, averageDayofInfectionTracingMatrix,
    xlabel="Test Sensitivity (1 - False Negative Rate)",
    ylabel="Tracing Efficiency",
    title="Average Day of Infection (Traced Cases and waiting time $(WaitBeforeTestTaken + WaitBeforeTestResult) days)",
    color=:viridis)

# Display both plots
plot(p1, p2, layout=(1,2), size=(1400,800))



# Create x-axis values (false negative test rates)
false_negative_rates = 0:0.02:1.0;

# Create the plot
plot(false_negative_rates, empiricalMean, 
    label="Empirical Mean", 
    linewidth=2, 
    xlabel="False Negative Test Rate", 
    ylabel="Effective Reproduction Number",
    xlim=(0,1),
    ylim=(0,R0),
    title="Wait Test: $(WaitBeforeTestTaken) days, Wait Test Result: $(WaitBeforeTestResult) days"
    )

plot!(false_negative_rates, theoreticalMeanArray, 
    label="Theoretical Mean", 
    xlim=(0,1),
    ylim=(0,R0),
    linewidth=2)

plot(false_negative_rates, empiricalMean, 
    label="Empirical Mean", 
    linewidth=2, 
    xlabel="False Negative Test Rate", 
    ylabel="Effective Reproduction Number",
    xlim=(0,1),
    )

plot!(false_negative_rates, theoreticalMeanArray, 
    label="Theoretical Mean", 
    xlim=(0,1),
    linewidth=2)

# Display the plot
display(plot!())

plot(false_negative_rates, (empiricalMean-theoreticalMeanArray), 
    linewidth=2, 
    label="Difference (Empirical - Theoretical)",
    xlabel="False Negative Test Rate", 
    ylabel="difference")


#test loops 
N = 100000
NumberOfChildrenToDraw = 0
I1 = 0.0
I2 = 0.0
global sumNumberOfChildrenToDraw = 0.0
sumofInfected_on_last_day = 0

N = MaximumAllowedInfected
τ = WaitBeforeTestTaken + WaitBeforeTestResult 
I1 = 0.0
I2 = 0.0
I3 = 0.0
AsymptomaticFractionOfInfected = 1
for test in 1:N    
    global GoalOfCountDown_individual = 2*drawLognormallyDistributedInteger();
    t_half = GoalOfCountDown_individual/2 + 1

    InfectiousnessDistribution = getInfectiousnessDistribution(GoalOfCountDown_individual, InfectiousProfile);
    L = length(InfectiousnessDistribution);

    # Convert to integer indices
    b1 = floor(Int, t_half);
    b2 = floor(Int, t_half + τ);

    if b2 <= L
        I1 += sum(InfectiousnessDistribution[b2+1:end]);
        I2 += sum(InfectiousnessDistribution[b1+1:b2]);
    else
        I1 += 0;
        I2 += sum(InfectiousnessDistribution[b1+1:end]);
    end
    #I3 += sum(InfectiousnessDistribution[1:b1]);
                 
    # Find out how many people this infectious node will infect.
    global R0OfNode = R0 * GoalOfCountDown_individual / (2*MeanOfLognormal);
    global NumberOfChildrenToDraw = drawNumberOfChildren(R0OfNode,OffspringDistribution)
    global WhenInfectedWillInfectOthers_individual = drawTimesWhenInfectedWillInfectOthers(NumberOfChildrenToDraw, GoalOfCountDown_individual, InfectiousProfile);
    global sumNumberOfChildrenToDraw += NumberOfChildrenToDraw
    if NumberOfChildrenToDraw != 0;
        sumofInfected_on_last_day += sum(GoalOfCountDown_individual .== WhenInfectedWillInfectOthers_individual)
    end
end
print("R_eff_theoretical: \t", R0 - R0*(1-AsymptomaticFractionOfInfected)*((1-ProbabilityFalseNegativeTest)*I1/N+I2/N), "\n")

println("Average number of children drawn:\t", sumNumberOfChildrenToDraw/N)


# Check distribution of infection times
observations = [];
R0 = 2;
infectlength = 10
for i in 1:100000
    k = drawNumberOfChildren(2,"poisson");
    push!(observations, drawTimesWhenInfectedWillInfectOthers(k,infectlength,"empirical")...);
end
histogram(observations, bins=0.5:1:infectlength+0.5, xlabel="Day of infection", ylabel="Frequency", title="Histogram of Infection Times", normalized=true)
scatter!(1:infectlength, getInfectiousnessDistribution(infectlength,"empirical"), 
         label="Infectiousness Distribution", alpha=0.5, 
         legend=:bottom)


#Test mean of t_infect
R0 = 10;
totalInfectionTime = 0.0;
numSamples = 100000;
totalSum = 0;
samples = [];
for i in 1:numSamples
    GoalOfCountDown_individual = 2*drawLognormallyDistributedInteger();
    R0OfNode = R0 * GoalOfCountDown_individual / (2*MeanOfLognormal);
    NumberOfChildrenToDraw = drawNumberOfChildren(R0OfNode,"poisson")
    if NumberOfChildrenToDraw == 0
        continue
    end
    if GoalOfCountDown_individual > 119
        continue
    end
    WhenInfectedWillInfectOthers_individual = drawTimesWhenInfectedWillInfectOthers(NumberOfChildrenToDraw,GoalOfCountDown_individual,"empirical")
    normalisedWhenWillInfect = WhenInfectedWillInfectOthers_individual ./ GoalOfCountDown_individual
    
    samples = push!(samples, normalisedWhenWillInfect...)

end

using Plots
histogram(samples, bins=20, xlabel="Normalized Infection Time", ylabel="Frequency", title="Histogram of Normalized Infection Times", alpha=0.7, normalize=:pdf)
scatter!(0:0.05:1, getInfectiousnessDistribution(20,"empirical"), 
         label="Infectiousness Distribution", alpha=0.5, 
         legend=:topright)

p_false = 0.1 #round(rand(), digits=2)
GoalOfCountDown_individual = 8 #2*drawLognormallyDistributedInteger();

weighted_p_test = getFalseNegativeProbabilityDistribution(p_false, getInfectiousnessDistribution(GoalOfCountDown_individual,"empirical"), GoalOfCountDown_individual)
print("mean is " ,mean(weighted_p_test))   
plot((1-p_false).*getInfectiousnessDistribution(GoalOfCountDown_individual,"empirical")*GoalOfCountDown_individual, xlabel="Weighted False Negative Probability", alpha=0.7, label = "Weighted test sensitivity profile", title = "p_false=$(p_false)  ")
plot!(weighted_p_test, xlabel="Adjusted Weighted test sensitivity Probability", alpha=0.7, label = "adjusted Weighted test sensitivity profile")
plot!(getInfectiousnessDistribution(GoalOfCountDown_individual,"empirical"), xlabel="infectious profile", alpha=0.7, label ="Infectious profile")

