function AdvanceInfectedOneTimestep(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,NodeCanTestPositive,TraceNodesChildren,MaximumAllowedInfected,R0,MeanOfLognormal,WaitBeforeTestResult,NumberOfRecovered,ProbabilityFalseNegativeTest,OffspringDistribution,InfectiousProfile)
    # This function advances all infectious nodes 1 step in their course of disease.

    # Inputs:
    # ---------
    #   StateOfNodes:                   Vector. Entry i shows state of node. 0 if Susceptible, 1 if Exposed, 2 if Infectious, 3 if Removed.
    #   CountUpToStateChange:           Vector. Entry i counts up towhen state of node i changes; Type: Vector
    #   GoalOfCountDown:                Vector. Entry i determines when state of node i changes. Type: Vector.
    #   WhenInfectedWillInfectOthers:   Vector of lists. List on entry i contains times after infection when node i will infect others.
    #   TestArrivalTimeOfNodes:         Vector. Each entry is a countdown to when a node will get tested. 0 at day when test will be taken.
    #   ResultArrivalTimeOfNodes:       Vector. Each entry is a countdown to arrival of test result. 0 at day of result arrival.
    #   NodeCanTestPositive:            Vector. Entry i is 1 if node i was infectious when it was tested. Otherwise 0.
    #   TraceNodesChildren              Vector. Entry i is 1 if node i's children should be traced.

    #   MaximumAllowedInfected:         Integer. Maximum number of nodes whose course of disease we follow.
    #   R0:                             Float. Basic reproduction number.
    #   MeanOfLognormal:                Float. Mean of lognormal distribution used for incubation period.
    #   WaitBeforeTestResult:           Integer. Number of days before a test result arrives after test is taken.
    #   NumberOfRecovered:              Integer. Number of nodes that have recovered from disease since start of simulation.
    #   ProbabilityFalseNegativeTest:   Float. Probability that a test of an infectious node will turn out to be negative.
    #   OffspringDistribution:          String. Specifies the offspring distribution of the simulation. 
    #   InfectiousProfile:              String. Specifies temporal profile of infectiousness.


    # Outputs:
    # ---------

    # Variable to check if all nodes are susceptible or removed.
    FoundNoInfectiousOrExposedNode = true;

    # Loop over nodes that might be infectious.
    for InfectedNode=1:MaximumAllowedInfected
        # Check if node is infectious or exposed.
        if StateOfNodes[InfectedNode] == 1 || StateOfNodes[InfectedNode] == 2
            # If it is, remember.
            FoundNoInfectiousOrExposedNode = false;

            # Advance counter to state change 1 step
            CountUpToStateChange[InfectedNode] += 1;

            # Check if node changes state.
            if CountUpToStateChange[InfectedNode] == GoalOfCountDown[InfectedNode]
                # Change state of node.

                
                if StateOfNodes[InfectedNode]==1
                    # If node is Exposed.

                    # Make node Infectious
                    StateOfNodes[InfectedNode],CountUpToStateChange[InfectedNode],GoalOfCountDown[InfectedNode],WhenInfectedWillInfectOthers[InfectedNode] = makeNodeInfectious(StateOfNodes[InfectedNode],CountUpToStateChange[InfectedNode],GoalOfCountDown[InfectedNode],WhenInfectedWillInfectOthers[InfectedNode],R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile);

                elseif StateOfNodes[InfectedNode]==2
                    # If node is Infectious
                    
                    # Make node Removed
                    StateOfNodes[InfectedNode],CountUpToStateChange[InfectedNode] = makeNodeRemoved(StateOfNodes[InfectedNode],CountUpToStateChange[InfectedNode]);

                    # Remember that 1 node recovered.
                    NumberOfRecovered +=1;

                end


            end
        end
        # Check if node is waiting for a test .
        StateOfNodes[InfectedNode],TestArrivalTimeOfNodes[InfectedNode],ResultArrivalTimeOfNodes[InfectedNode],NodeCanTestPositive[InfectedNode] = AdvanceTestWaitOneStep(StateOfNodes[InfectedNode],TestArrivalTimeOfNodes[InfectedNode],ResultArrivalTimeOfNodes[InfectedNode],NodeCanTestPositive[InfectedNode],WaitBeforeTestResult)

        # Check if node is waiting for a test result.
        StateOfNodes[InfectedNode],CountUpToStateChange[InfectedNode],GoalOfCountDown[InfectedNode],TestArrivalTimeOfNodes[InfectedNode],ResultArrivalTimeOfNodes[InfectedNode],NodeCanTestPositive[InfectedNode],TraceNodesChildren[InfectedNode],NumberOfRecovered = AdvanceResultWaitOneStep(StateOfNodes[InfectedNode],CountUpToStateChange[InfectedNode],GoalOfCountDown[InfectedNode],TestArrivalTimeOfNodes[InfectedNode],ResultArrivalTimeOfNodes[InfectedNode],NodeCanTestPositive[InfectedNode],TraceNodesChildren[InfectedNode],NumberOfRecovered,ProbabilityFalseNegativeTest)

    end
    return StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,NodeCanTestPositive,TraceNodesChildren,WaitBeforeTestResult,NumberOfRecovered,FoundNoInfectiousOrExposedNode
end

function AdvanceResultWaitOneStep(StateOfNodes_specific,CountUpToStateChange_specific,GoalOfCountDown_specific,TestArrivalTimeOfNodes_specific,ResultArrivalTimeOfNodes_specific,NodeCanTestPositive_specific,TraceNodesChildren_specific,NumberOfRecovered,ProbabilityFalseNegativeTest)

    # Check if node is waiting for a test result
    if ResultArrivalTimeOfNodes_specific >=0
        # If node is waiting for a test result, check if result arrives today.
        if ResultArrivalTimeOfNodes_specific == 0 
            # If result arrives today draw a random number to check if test is a false negative.
            ResultRandomVariable = rand()+0;

            # If result arrives today, check whether result is false negative and whether node can even test positive.
            if NodeCanTestPositive_specific==1 && ResultRandomVariable>ProbabilityFalseNegativeTest
                # If node has not already recovered, add one node to recovered population
                if StateOfNodes_specific != 3
                    NumberOfRecovered +=1;
                end
                # Node tested positive. Trace its children.
                TraceNodesChildren_specific =1+0;

                # Node gets isolated (in model this is done by letting node recover).
                StateOfNodes_specific = 3+0;

                # Ignore any counting down to state change, tests or test results.
                CountUpToStateChange_specific = -1 +0;
                GoalOfCountDown_specific = -9 +0;
                TestArrivalTimeOfNodes_specific = -9+0;
                ResultArrivalTimeOfNodes_specific = -9+0;


            else 
                # Node did not test positive. Register that result has arrived.
                ResultArrivalTimeOfNodes_specific = -9+0;

            end
            # Result was delivered. Reset whether node can test positive.
            NodeCanTestPositive_specific = 0;

        else 
            # If Result was not delivered this time step, advance time 1 step.
            ResultArrivalTimeOfNodes_specific -=1;
        end


    end

    return StateOfNodes_specific,CountUpToStateChange_specific,GoalOfCountDown_specific,TestArrivalTimeOfNodes_specific,ResultArrivalTimeOfNodes_specific,NodeCanTestPositive_specific,TraceNodesChildren_specific,NumberOfRecovered
end

function AdvanceTestWaitOneStep(StateOfNodes_specific,TestArrivalTimeOfNodes_specific,ResultArrivalTimeOfNodes_specific,NodeCanTestPositive_specific,WaitBeforeTestResult)
    # This function checks if a node is waiting for having a test taken. If test is due, it starts the node's wait for a test result.

    # Inputs
    # --------
    # StateOfNodes_specific:                State of the node in question. 0 If susceptible, 1 if exposed, 2 if infectious, 3 if recovered.
    # TestArrivalTimeOfNodes_specific:      Countdown to when a node will get tested. 0 at day when test will be taken.
    # ResultArrivalTimeOfNodes_specific:    Countdown to arrival of test result. 0 at day of result arrival.
    # NodeCanTestPositive_specific:         Integer. 1 if node was infectious when it was tested. Otherwise 0.
    # WaitBeforeTestResult:                 Integer. Number of days before a test result arrives after test is taken.

    # Outputs
    # --------
    # StateOfNodes_specific:                Integer. State of the node in question. 0 If susceptible, 1 if exposed, 2 if infectious, 3 if recovered.
    # TestArrivalTimeOfNodes_specific:      Integer. Countdown to when a node will get tested. 0 at day when test will be taken.
    # ResultArrivalTimeOfNodes_specific:    Integer. Countdown to arrival of test result. 0 at day of result arrival.
    # NodeCanTestPositive_specific:         Integer. 1 if node was infectious when it was tested. Otherwise 0.


    # Check if node is waiting for a test
    if TestArrivalTimeOfNodes_specific >= 0

        # If waiting, check if test is taken today.
        if TestArrivalTimeOfNodes_specific == 0
            # If test is taken today, start node's wait for a test result.
            ResultArrivalTimeOfNodes_specific = WaitBeforeTestResult+0;

            # Check whether this test can be positive. Only can if node is infectious at time of test.
            if StateOfNodes_specific == 2
                NodeCanTestPositive_specific = 1 +0;
            end
        end
        # Advance node's test count down one.
        TestArrivalTimeOfNodes_specific -=1;
    end

    return StateOfNodes_specific,TestArrivalTimeOfNodes_specific,ResultArrivalTimeOfNodes_specific,NodeCanTestPositive_specific
    
end

function AppendLineToFile(FileDestination,StringToSave)
    # This function saves a string as a new line in a file.
    # Inputs:
    # - FileDestination:    String. Path to where line should be saved.
    # - StringToSave:       String. The string that should be saved to the file.
    open(FileDestination,"a") do File
    
        write(File,string("\n",StringToSave))
    end

end

function drawGeometricInteger(R0OfNode)
    # Return an integer drawn from a Geometric distribution with mean R0OfNode
    # Inputs 
    # --------
    # R0OfNode:     Float. Mean of Geometric distribution

    # Outputs
    # --------
    # DrawnInteger  Integer. Drawn from Geometric distribution

    RandomFloat = rand()+0;
    DrawnInteger = -1;

    CumulativeGeometricDistribution =0 ;
    DrawnIntegerFound = false;
    while DrawnIntegerFound==false
        DrawnInteger += 1;

        CumulativeGeometricDistribution += evaluateGeometricDistribution(R0OfNode,DrawnInteger)

        if CumulativeGeometricDistribution >= RandomFloat
            DrawnIntegerFound = true
            break
        end
    end
    return DrawnInteger
end



function drawLognormallyDistributedInteger()
    # Output:
    # An integer drawn from a lognormal distribution.


    # Retrieve lognormal distribution
	LognormalDistribution = getLognormalDistribution();
	RandomFloat = rand()+0;

	SumOfLognormallyDistributedValues =0;

	Result = length(LognormalDistribution) + 0;
	for Day = 1:length(LognormalDistribution)-1

		SumOfLognormallyDistributedValues += LognormalDistribution[Day];
		if SumOfLognormallyDistributedValues >= RandomFloat
			Result = Day;
			break

        end

	end

	return Result


end


function drawNumberOfChildren(R0OfNode,OffspringDistribution)

    if (OffspringDistribution=="poisson")

        NumberOfChildren = drawPoissonInteger(R0OfNode)

    elseif (OffspringDistribution=="geometric")
        NumberOfChildren = drawGeometricInteger(R0OfNode)

    end
    return NumberOfChildren
end

function drawPoissonInteger(R0OfNode)
    # Return an integer drawn from a Poisson distribution with mean R0OfNode
    # Inputs 
    # --------
    # R0OfNode:     Float. Mean of Poisson distribution

    # Outputs
    # --------
    # DrawnInteger  Integer. Drawn from Poisson distribution

    RandomFloat = rand()+0;
    DrawnInteger = -1;

    CumulativePoissonDistribution =0 ;
    DrawnIntegerFound = false;
    while DrawnIntegerFound==false
        DrawnInteger += 1;

        CumulativePoissonDistribution += evaluatePoissonDistribution(R0OfNode,DrawnInteger)

        if CumulativePoissonDistribution >= RandomFloat
            DrawnIntegerFound = true
            break
        end
    end
    return DrawnInteger
end

function drawTimesWhenInfectedWillInfectOthers(NumberOfChildrenToDraw,CounterGoalOfNode,InfectiousProfile)
    # Creates an array containing at what times an infectious node will infect others.
    # Inputs:
    # ---------
    # - NumberOfChildrenToDraw:     Integer. The number of people the node will infect.
    # - CounterGoalOfNode:          Integer. The number of days the node will be infectious.
    # - InfectiousProfile:          String. Name of the profile of infectiousness.

    # Outputs:
    # ---------
    # TimesOfInfection              Array of integers. Contains times at which the node infects other people.

    # Make array to contain infection times.
    TimesOfInfection = [];

    # Get infectiousness profile of the period that node is infectious.
    InfectiousnessDistribution = getInfectiousnessDistribution(CounterGoalOfNode,InfectiousProfile)

    for ChildNumber =1:NumberOfChildrenToDraw
        if InfectiousProfile=="empirical" || InfectiousProfile=="FlatSkewed"
            RandomFloat = rand()+0;
            Day =0;

            FoundDayOfInfection = false;
            SumInfectiousProfile = 0;
            while FoundDayOfInfection == false
                Day +=1;
                SumInfectiousProfile += InfectiousnessDistribution[Int(Day)];

                if SumInfectiousProfile>= RandomFloat
                    FoundDayOfInfection = true
                    break
                end

            end
            TimesOfInfection = push!(TimesOfInfection,Day);

        else
            print("USING INFECTIOUS PROFILE THAT WAS NOT IMPLEMENTED")
        end
    end
    return TimesOfInfection
end

function evaluateGeometricDistribution(R0OfNode,DrawnInteger)
    return (R0OfNode/(1+R0OfNode))^DrawnInteger / (R0OfNode + 1)
end

function evaluatePoissonDistribution(R0OfNode,DrawnInteger)
    if DrawnInteger > 20
        factorial_computed = factorial(big(DrawnInteger));
    else
        factorial_computed = factorial(DrawnInteger);
    end
    return R0OfNode^DrawnInteger/factorial_computed*exp(-R0OfNode)
end
function getCorrectedEmpiricalInfectiousness() 

	CorrectedEmpiricalInfectiousness = [0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000006, 0.000000000400, 0.000000014828, 0.000000327126, 0.000004598193, 0.000043551427, 0.000291074036, 0.001426631003, 0.005296911006, 0.015315428062, 0.035311746313, 0.066265116634, 0.103030125971, 0.134813059932, 0.150505975534, 0.145114140518, 0.122150668664, 0.090637282629, 0.059800573129, 0.035357288295, 0.018866386876, 0.009143443810, 0.004048236058, 0.001646108958, 0.000617722128, 0.000214881890, 0.000069574929, 0.000021046826, 0.000005969123, 0.000001592286, 0.000000400693, 0.000000095386, 0.000000021536, 0.000000004623, 0.000000000945, 0.000000000185, 0.000000000034, 0.000000000006, 0.000000000001, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000];
	return CorrectedEmpiricalInfectiousness # entry 29 is symptom onset.
end

function getFlatSkewedInfectiousness()
    FlatSkewedInfectiousness = [0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.011204481792717087,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435,0.0056022408963585435];
    return FlatSkewedInfectiousness
end
function getInfectiousnessDistribution(CounterGoalOfNode,InfectiousProfile)

	InfectiousnessDistribution = zeros(Int(CounterGoalOfNode))#fill(Int[], Int(CounterGoalOfNode),1);

    if (InfectiousProfile=="empirical")
		EmpiricalInfectiousness = getCorrectedEmpiricalInfectiousness();
    elseif (InfectiousProfile=="FlatSkewed")
        EmpiricalInfectiousness = getFlatSkewedInfectiousness();
    else
        print("USING INFECTIOUS PROFILE THAT WAS NOT IMPLEMENTED")
    end

    DenominatorSum = 0;
    OffsetDays = Int(max(floor((length(EmpiricalInfectiousness)-1)/2)-floor(CounterGoalOfNode/2), 0));

    for Day = 1:Int(CounterGoalOfNode)
        InfectiousnessDistribution[Day] += EmpiricalInfectiousness[OffsetDays+Day]
        DenominatorSum += EmpiricalInfectiousness[OffsetDays+Day];

    end

    for Day = 1:Int(CounterGoalOfNode)
        InfectiousnessDistribution[Day] /= DenominatorSum;
    end        



    return InfectiousnessDistribution
end


function getInitialConditionsOfSimulation(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,InitialNumberOfInfected,R0,MeanOfLognormalDistribution,OffspringDistribution,InfectiousProfile)
    # Initialized 4 of state vectors.
    # Inputs:
    # ----------
    #   StateOfNodes:                   Vector with state of node i on entry i; type: Vector
    #   CountUpToStateChange:           Vector. Entry i counts up towhen state of node i changes; Type: Vector
    #   GoalOfCountDown:                Vector. Entry i determines when state of node i changes. Type: Vector.
    #   WhenInfectedWillInfectOthers:   Vector of lists. List on entry i contains times after infection when node i will infect others.
    #   InitialNumberOfInfected:        Integer. How many nodes are infected at the beginning of a simulation.
    #   R0:                             Float. Basic reproduction number.
    #   MeanOfLognormal:                Float. Mean of lognormal distribution used for incubation period.
    #   OffspringDistribution:          String. Specifies the offspring distribution of the simulation. 
    #   InfectiousProfile:              String. Specifies temporal profile of infectiousness.

    # Outputs:
    # ----------
    #   StateOfNodes:           Vector with state of node i on entry i; type: Vector
    #   CountUpToStateChange:   Vector. Entry i counts up towhen state of node i changes; Type: Vector
    #   GoalOfCountDown:        Vector. Entry i determines when state of node i changes. Type: Vector.
    #   WhenInfectedWillInfectOthers:   Vector of lists. List on entry i contains times after infection when node i will infect others.

    # Seeds start out Infectious.
    StateOfNodes[1:InitialNumberOfInfected] .+= 2;

    # Seeds start counting up to state change.
    CountUpToStateChange[1:InitialNumberOfInfected] .= 0;

    for SeedNumber = 1:InitialNumberOfInfected
        # Draw goal for count of all seeds    
        GoalOfCountDown[SeedNumber] += 2*drawLognormallyDistributedInteger();

        # We expect people with longer infectious periods to infect more...
        R0OfNode = R0 * GoalOfCountDown[SeedNumber] / (2*MeanOfLognormalDistribution) ; # Assume infectivity is proporsional to length of infectious period. Factor of 2 because Goal is 2 times Lognormally drawn integer..


        # Draw times when seeds will infect others.
        NumberOfChildrenThisNodeInfects = drawNumberOfChildren(R0OfNode,OffspringDistribution);
        WhenInfectedWillInfectOthers[SeedNumber] = drawTimesWhenInfectedWillInfectOthers(NumberOfChildrenThisNodeInfects, GoalOfCountDown[SeedNumber], InfectiousProfile);

    end

    return StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers
end

function getLognormalDistribution()
	# Incubation time lognormal distribution
	LognormalDistribution = [0.003928899789, 0.068966842356, 0.151229922834, 0.177164768798, 0.159473206432, 0.126225212682, 0.093377539613, 0.066601193262, 0.046592864872, 0.032292432894, 0.022307832766, 0.015417789810, 0.010686160309, 0.007438793831, 0.005205546862, 0.003663990068, 0.002594779101, 0.001849127194, 0.001326071478, 0.000956932785, 0.000694814444, 0.000507543367, 0.000372931527, 0.000275591462, 0.000204790215, 0.000152997488, 0.000114898935, 0.000086722192, 0.000065774043, 0.000050120846, 0.000038366621, 0.000029498110, 0.000022775974, 0.000017657978, 0.000013744467, 0.000010739444, 0.000008422656, 0.000006629465, 0.000005236258, 0.000004149825, 0.000003299588, 0.000002631891, 0.000002105782, 0.000001689883, 0.000001360066, 0.000001097709, 0.000000888391, 0.000000720903, 0.000000586509, 0.000000478372, 0.000000391131, 0.000000320564, 0.000000263340, 0.000000216822, 0.000000178915, 0.000000147954, 0.000000122607, 0.000000101811, 0.000000084712];

	return LognormalDistribution
end

function getMeanOfLognormalDistribution()
	MeanOfLognormalDistribution =0;
	LognormalDistribution = getLognormalDistribution();
	for Day = 1:length(LognormalDistribution)
		MeanOfLognormalDistribution += (Day) * LognormalDistribution[Day];
    end
	return MeanOfLognormalDistribution
end

function InfectNodesOnThisTimestep(StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,MaximumAllowedInfected,NumberOfInfected,R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile)

    # Inputs
    # --------
    #   StateOfNodes:                   Vector. Entry i shows state of node. 0 if Susceptible, 1 if Exposed, 2 if Infectious, 3 if Removed.
    #   CountUpToStateChange:           Vector. Entry i counts up towhen state of node i changes; Type: Vector
    #   GoalOfCountDown:                Vector. Entry i determines when state of node i changes. Type: Vector.
    #   WhenInfectedWillInfectOthers:   Vector of lists. List on entry i contains times after infection when node i will infect others.
    #   ListOfChildren:
    #   TestArrivalTimeOfNodes:         Vector. Each entry is a countdown to when a node will get tested. 0 at day when test will be taken.
    #   ResultArrivalTimeOfNodes:       Vector. Each entry is a countdown to arrival of test result. 0 at day of result arrival.
    
    #   MaximumAllowedInfected:         Integer. Maximum number of nodes whose course of disease we follow.
    #   NumberOfInfected   
    #   R0:                             Float. Basic reproduction number (mean).
    #   MeanOfLognormal:                Float. Mean of lognormal distribution defining incubation time.
    #   OffspringDistribution:          String. Specifies the offspring distribution of the simulation. 
    #   InfectiousProfile:              String. Specifies temporal profile of infectiousness.



    # Outputs
    # --------
    #   StateOfNodes:                   Vector. Entry i shows state of node. 0 if Susceptible, 1 if Exposed, 2 if Infectious, 3 if Removed.
    #   CountUpToStateChange:           Vector. Entry i counts up towhen state of node i changes; Type: Vector
    #   GoalOfCountDown:                Vector. Entry i determines when state of node i changes. Type: Vector.
    #   WhenInfectedWillInfectOthers:   Vector of lists. List on entry i contains times after infection when node i will infect others.
    #   ListOfChildren:
    #   NumberOfInfected

    # Loop over nodes that might be infected.
    for InfectedNode=1:MaximumAllowedInfected
        # Check if node is 1) Infectious, 2) Waiting for a test, 3) Waiting for a test result
        if StateOfNodes[InfectedNode]==2 && TestArrivalTimeOfNodes[InfectedNode]<0 && ResultArrivalTimeOfNodes[InfectedNode]<0
            # Loop over node's coming children
            for NodesChildNumber=1:length(WhenInfectedWillInfectOthers[InfectedNode])
                # Check if this child is due today.
                if WhenInfectedWillInfectOthers[InfectedNode][NodesChildNumber] == CountUpToStateChange[InfectedNode]
                    # Nodes gets infected.
                    NumberOfInfected +=1;
                    # Check if more infected nodes are allowed.
                    if NumberOfInfected <= MaximumAllowedInfected
                        # BEWARE: If exposed state is introduced, the following will make trouble.
                        # Make the node infectious
                        StateOfNodes[NumberOfInfected],CountUpToStateChange[NumberOfInfected],GoalOfCountDown[NumberOfInfected],WhenInfectedWillInfectOthers[NumberOfInfected] = makeNodeInfectious(StateOfNodes[NumberOfInfected],CountUpToStateChange[NumberOfInfected],GoalOfCountDown[NumberOfInfected],WhenInfectedWillInfectOthers[NumberOfInfected],R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile);
                    end
                    # Add newly infected to infecting node's list of children. Used for tracing. 
                    # (have to copy in order to not add the node to _all_ lists in the list)
                    ListOfChildren_specific = copy(ListOfChildren[InfectedNode]);
                    ListOfChildren_specific = push!(ListOfChildren_specific,NumberOfInfected+0);
                    ListOfChildren[InfectedNode] = copy(ListOfChildren_specific);
                end
            end
        end
    end

    return StateOfNodes,CountUpToStateChange,GoalOfCountDown,WhenInfectedWillInfectOthers,ListOfChildren,NumberOfInfected
end

function makeNodeInfectious(StateOfNodes_individual,CountUpToStateChange_individual,GoalOfCountDown_individual,WhenInfectedWillInfectOthers,R0,MeanOfLognormal,OffspringDistribution,InfectiousProfile)
    # This function changes the state of a node from Exposed to Infectious.

    # Inputs
    # --------
    # StateOfNodes_individual:                   Integer. State of node in question. Is node's entry in StateOfNodes vector.
    # CountUpToStateChange_individual:          Integer. Counter up to next time the node changes state. Is node's entry in CountUpToStateChange vector.
    # GoalOfCountDown_individual:               Integer. When the node will change state next time.
    # WhenInfectedWillInfectOthers_individual:  List. Lists times after infection that node will infect other nodes.
    # R0:                                       Float. Basic reproduction number (mean).
    # MeanOfLognormal:                          Float. Mean of lognormal distribution defining incubation time.
    # OffspringDistribution:                    String. Specifies the offspring distribution of the simulation. 
    # InfectiousProfile:                        String. Specifies temporal profile of infectiousness.

    # Outputs
    # ---------

    # First, make state to 2:
    StateOfNodes_individual =2+0;

    # Reset count down to next change of node state.
    CountUpToStateChange_individual =0;

    # Draw time for node's next state change.
    GoalOfCountDown_individual = 2*drawLognormallyDistributedInteger(); ;

    # Find out how many people this infectious node will infect.
    R0OfNode = R0 * GoalOfCountDown_individual / (2*MeanOfLognormal);
    NumberOfChildrenToDraw = drawNumberOfChildren(R0OfNode,OffspringDistribution);
    WhenInfectedWillInfectOthers_individual = drawTimesWhenInfectedWillInfectOthers(NumberOfChildrenToDraw, GoalOfCountDown_individual, InfectiousProfile);

    return StateOfNodes_individual,CountUpToStateChange_individual,GoalOfCountDown_individual,WhenInfectedWillInfectOthers_individual

end


function makeNodeRemoved(StateOfNodes_individual,CountUpToStateChange_individual)
    # This function changes the state of a node from Infectious to Removed.

    # Inputs
    # --------
    # StateOfNodes_individual:          Integer. State of node in question. Is node's entry in StateOfNodes vector.
    # CountUpToStateChange_individual:          Integer. Counter up to next time the node changes state. Is node's entry in CountUpToStateChange vector.

    # Outputs
    # ---------

    # Advance state 1:
    StateOfNodes_individual +=1;

    # Stop counter to count up to another change of states.
    CountUpToStateChange_individual = -1;

    return StateOfNodes_individual,CountUpToStateChange_individual

end


function TraceNode(StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren,Asymptomatic,MaximumAllowedInfected,WaitBeforeTestTaken,ProbabilityChildIsTraced)

    # This function traces children of nodes that tested positive this time step. Also orders test for nodes that got symptomatic this time step.

    # Inputs
    # --------
    #   StateOfNodes:                   Vector. Entry i shows state of node. 0 if Susceptible, 1 if Exposed, 2 if Infectious, 3 if Removed.
    #   CountUpToStateChange:           Vector. Entry i counts up towhen state of node i changes; Type: Vector
    #   GoalOfCountDown:                Vector. Entry i determines when state of node i changes. Type: Vector.
    #   ListOfChildren:                 List of Lists. List i contains IDs of nodes that node i infected.
    #   TestArrivalTimeOfNodes:         Vector. Each entry is a countdown to when a node will get tested. 0 at day when test will be taken.
    #   ResultArrivalTimeOfNodes:       Vector. Each entry is a countdown to arrival of test result. 0 at day of result arrival.
    #   TraceNodesChildren:             Vector: Entry i is 1 if node i is being contact traced.
    #   Asymptomatic:                   Vector. Entry i is 1 if node i will not develop symptoms.

    #   MaximumAllowedInfected:         Integer. Number of nodes that program follows through course of disease.
    #   WaitBeforeTestTaken:            Integer. Number of days between a test is ordered and it is taken.
    #   ProbabilityChildIsTraced:       Float. Probability that contact tracing will successfully identify each child of a node being contact traced.

    # Outputs
    # ---------


    for FocalNode=1:MaximumAllowedInfected
        # Check if nodes gets symptomatic
        if StateOfNodes[FocalNode] == 2 && CountUpToStateChange[FocalNode]>0 && CountUpToStateChange[FocalNode] == floor(GoalOfCountDown[FocalNode]/2+1) && TestArrivalTimeOfNodes[FocalNode]<0 && ResultArrivalTimeOfNodes[FocalNode]<0 && Asymptomatic[FocalNode]==0
            # If it does and is not waiting for test or result, order test.
            TestArrivalTimeOfNodes[FocalNode]=WaitBeforeTestTaken+0;
        end

        # Check whether Node's children are to be traced this time step.
        if TraceNodesChildren[FocalNode]==1 
            for NodesChildNumber=1:length(ListOfChildren[FocalNode])
                # Draw random number. This decides if child is successfully traced.
                TraceRandomVariable = rand()+0;
                if TraceRandomVariable < ProbabilityChildIsTraced
                    # Child was traced. If child is one of the nodes we are following and it is not already being tested, order a test.
                    IDOfTracedChild = ListOfChildren[FocalNode][NodesChildNumber] + 0;
                    if IDOfTracedChild<=MaximumAllowedInfected
                        if TestArrivalTimeOfNodes[IDOfTracedChild]<0 && ResultArrivalTimeOfNodes[IDOfTracedChild]<0
                            TestArrivalTimeOfNodes[IDOfTracedChild] = WaitBeforeTestTaken + 0;
                        end
                    end
                
                end

            end

            TraceNodesChildren[FocalNode] = 0+0;
        end

    end
    return StateOfNodes,CountUpToStateChange,GoalOfCountDown,ListOfChildren,TestArrivalTimeOfNodes,ResultArrivalTimeOfNodes,TraceNodesChildren
end
