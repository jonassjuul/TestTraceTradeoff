using DelimitedFiles
include("Functions_main.jl")


# Loop through wait times 0-5
matrices = []
filename = "code/OutputsMorten/averageDayofInfectionMatrix_wait0.txt"
referencematrix = zeros(10,10)

for wait in 0:5
    filename = "code/OutputsMorten/averageDayofInfectionTracingMatrix_wait$(wait).txt"
    referencefilename = "code/OutputsMorten/averageDayofInfectionMatrix_wait$(wait).txt"
    # Check if file exists before loading
    if isfile(filename)
        matrix = readdlm(filename)
        matrixref = readdlm(referencefilename)
        println("Loaded matrix for wait time $wait with size $(size(matrix))")
        # Store in a list - initialize outside the loop if needed

        matrix = matrix[2:end, 1:end-1]
        matrixref = matrixref[2:end, 1:end-1]/6

        push!(matrices, matrix)
        referencematrix += matrixref
        
    else
        println("File not found: $filename")
    end
end

#%% Plotting difference matrices

figurePath = "code/figuresMorten/"

matrices_diff = []
for (i, matrix) in enumerate(matrices)
    diff = matrix - referencematrix
    push!(matrices_diff, diff)
    println("Computed difference matrix for wait time $(i-1) with size $(size(diff))")
end

# Calculate global min and max for consistent colormap
global_min = minimum(minimum(diff) for diff in matrices_diff)
global_max = maximum(maximum(diff) for diff in matrices_diff)
false_negative_rates = 0.1:0.1:1.0;
trace_effeciency_rates = 0.1:0.1:1.0;
# Now plot the difference matrices
using Plots
plot_titles = ["Wait time = $(i-1) days" for i in 1:length(matrices_diff)]
p1 = plot(layout = (2, 3), size=(1200, 800), plot_title="Difference in Average Infection Day (Tracing - No Tracing)")
for (i, diff_matrix) in enumerate(matrices_diff)
    heatmap!(p1[i],false_negative_rates, trace_effeciency_rates, diff_matrix, title=plot_titles[i], xlabel="False Negative Rate", ylabel="Tracing Efficiency", 
             colorbar_title="Difference in Avg Infection Day", clims=(global_min, global_max))
end
display(p1)

savefig(p1, figurePath*"differenceInAvgInfectionday.png")
plot_titles = ["Wait time = $(i-1) days" for i in 1:length(matrices_diff)]
p4 = plot(layout = (3, 3), size=(1200, 800), plot_title="Infection Time Given Tracing")
global_min = minimum(minimum(matrix) for matrix in matrices)
global_max = maximum(maximum(matrix) for matrix in matrices)
for (i, matrix) in enumerate(matrices)
    heatmap!(p4[i],false_negative_rates, trace_effeciency_rates, matrix, title=plot_titles[i], xlabel="False Negative Rate", ylabel="Tracing Efficiency", 
             colorbar_title="Infection time given tracing with wait $(i-1) days", clims=(global_min, global_max))
end

# Hide empty subplots 7 and 8
plot!(p4[7], framestyle=:none, showaxis=false, grid=false, ticks=false)
plot!(p4[8], framestyle=:none, showaxis=false, grid=false, ticks=false)

heatmap!(p4[9],false_negative_rates, trace_effeciency_rates, referencematrix, title="reference matrix with no tracing", xlabel="False Negative Rate", ylabel="Tracing Efficiency", 
             colorbar_title="average infection time")
display(p4)

savefig(p4, figurePath*"infectionTimeGivenTracing.png")

matricesTracingNumber = []
matricesTracingProb = []
matricesPandemicEnded = []

for wait in 0:5
    filename = "code/OutputsMorten/averagenumbertracedMatrix_wait$(wait).txt"
    
    # Check if file exists before loading
    if isfile(filename)
        matrix = readdlm(filename)
        println("Loaded matrix for wait time $wait with size $(size(matrix))")
        push!(matricesTracingNumber, matrix)
    else
        println("File not found: $filename")
    end
    
    filename = "code/OutputsMorten/averageproportiontracedMatrix_wait$(wait).txt"
    
    # Check if file exists before loading
    if isfile(filename)
        matrix = readdlm(filename)
        println("Loaded matrix for wait time $wait with size $(size(matrix))")
        push!(matricesTracingProb, matrix)
    else
        println("File not found: $filename")
    end

    filename = "code/OutputsMorten/averagepandemicended_wait$(wait).txt"
    
    # Check if file exists before loading
    if isfile(filename)
        matrix = readdlm(filename)
        println("Loaded matrix for wait time $wait with size $(size(matrix))")
        push!(matricesPandemicEnded, matrix)
    else
        println("File not found: $filename")
    end
end

false_negative_rates = 0.0:0.1:1.0;
trace_effeciency_rates = 0.0:0.1:1.0;

# Calculate global min and max for each matrix type
global_min_num = minimum(minimum(matrix) for matrix in matricesTracingNumber)
global_max_num = maximum(maximum(matrix) for matrix in matricesTracingNumber)
global_min_prob = 0
global_max_prob = 1
global_min_pandemic = 0
global_max_pandemic = 1
Pasymp = 0.3
theoreticaltraceefficiencyMatrix = zeros(length(trace_effeciency_rates), length(false_negative_rates));
for (i, tracingefficiency) in enumerate(trace_effeciency_rates)
    for (j, falsenegativerate) in enumerate(false_negative_rates)
        theoreticaltraceefficiencyMatrix[i, j] = (1-Pasymp) * (1 - falsenegativerate) * tracingefficiency/(1-Pasymp * (1 - falsenegativerate) * tracingefficiency);

    end
end

heatmap(false_negative_rates, trace_effeciency_rates, theoreticaltraceefficiencyMatrix, 
        title="Theoretical Tracing Efficiency", xlabel="False Negative Rate", 
        ylabel="Tracing Efficiency", colorbar_title="Theoretical Trace Efficiency")
# Plot both matrices in pairs (6 rows, 3 columns for 6 wait times)
p2 = plot(layout = (6, 4), size=(1600, 1600))
for i in 1:length(matricesTracingNumber)
    wait_time = i - 1
    
    heatmap!(p2[i, 1], false_negative_rates, trace_effeciency_rates, matricesTracingNumber[i], 
             title="Number Traced (Wait: $(wait_time))", xlabel="False Negative Rate", 
             ylabel="Tracing Efficiency", colorbar_title="Number Traced", 
             clims=(global_min_num, global_max_num))
    
    # Second column: Tracing Probabilities
    heatmap!(p2[i, 2], false_negative_rates, trace_effeciency_rates, matricesTracingProb[i], 
             title="Proportion Traced empirical (Wait: $(wait_time))", xlabel="False Negative Rate", 
             ylabel="Tracing Efficiency", colorbar_title="Proportion Traced", 
             clims=(global_min_prob, global_max_prob))
    # Third column: Pandemic Ended
    heatmap!(p2[i, 3], false_negative_rates, trace_effeciency_rates, matricesPandemicEnded[i], 
             title="Pandemic Ended (Wait: $(wait_time))", xlabel="False Negative Rate", 
             ylabel="Tracing Efficiency", colorbar_title="Pandemic Ended", 
             clims=(global_min_pandemic, global_max_pandemic))
    heatmap!(p2[i, 4], false_negative_rates, trace_effeciency_rates, matricesTracingProb[i]-theoreticaltraceefficiencyMatrix, 
             title="Diff Empirical - Theoretical", xlabel="False Negative Rate", 
             ylabel="Tracing Efficiency")
end
display(p2)
plot!(p2, plot_title = "Total number of traced nodes for Different Wait Times (empirical results from simulations)")
savefig(p2, figurePath*"amountOfTracedPeopleForEachScenario.png")


Pasymp = 0.3
MeanOfLognormal = getMeanOfLognormalDistribution();

equation3Matrices = []
boolean_matrices = []
for wait in 0:5
    matrix = zeros(length(trace_effeciency_rates)-1, length(false_negative_rates)-1)
    for (i, tracingefficiency) in enumerate(trace_effeciency_rates[2:end])
        for (j, falsenegativerate) in enumerate(false_negative_rates[1:end-1])
            p = (1-Pasymp) * (1 - falsenegativerate) * tracingefficiency/(1-Pasymp * (1 - falsenegativerate) * tracingefficiency);
            matrix[i, j] = MeanOfLognormal + 1 - p/(1-p)*(matrices[wait+1][i, j]+wait)
            
        end
    end
    push!(equation3Matrices, matrix)
    push!(boolean_matrices, matrix .<= 0)
end
# Now plot theoretical trace efficiency heatmap
p5 = plot(layout = (6, 3), size=(1600, 1600))
globalmin = minimum(minimum(matrix) for matrix in equation3Matrices)
globalmax = maximum(maximum(matrix) for matrix in equation3Matrices)
for (i, wait) in enumerate(0:5)
    heatmap!(p5[i,1], false_negative_rates[1:end-1], trace_effeciency_rates[2:end], equation3Matrices[i], 
             title="Theoretical Infection Time (Wait: $(wait))", xlabel="False Negative Rate", 
             ylabel="Tracing Efficiency", colorbar_title="Theoretical Infection Time", clims=(globalmin, globalmax))
    heatmap!(p5[i,2], false_negative_rates[1:end-1], trace_effeciency_rates[2:end], boolean_matrices[i], 
             title="Infection Time =< 0 (Wait: $(wait))", xlabel="False Negative Rate", 
             ylabel="Tracing Efficiency", colorbar_title="Infection Time > 0")
    heatmap!(p5[i,3], false_negative_rates, trace_effeciency_rates, matricesPandemicEnded[i], 
             title="ended pandemics in simulations (Wait: $(wait))", xlabel="False Negative Rate", 
             ylabel="Tracing Efficiency", colorbar_title="proportion of ended pandemics")
end


display(p5)
plot!(p5, plot_title = "Evaluation of Equation 3 for Different Wait Times")
savefig(p5, figurePath*"evaluationOfEq3.png")


empiricalMeans = []
theoreticalMeans = []

for wait in 0:5
    filename_emp = "code/OutputsMorten/empiricalMean_wait$(wait)NumExp50.txt"
    filename_theo = "code/OutputsMorten/theoreticalMean_wait$(wait)NumExp50.txt"
    
    if isfile(filename_emp) && isfile(filename_theo)
        empiricalMean = readdlm(filename_emp)
        theoreticalMean = readdlm(filename_theo)
        push!(empiricalMeans, empiricalMean)
        push!(theoreticalMeans, theoreticalMean)
        println("Loaded empirical and theoretical means for wait time $wait")
    else
        println("Files not found for wait time $wait")
    end
end

false_negative_rates = 0.0:0.1:1.0
trace_effeciency_rates = 0.0:0.1:1.0

# Create plot with 6 rows (one for each wait time) and 3 columns
using Plots

eq3 = false  # Set to true to use Eq3 boolean mask, false for pandemic ended mask
p6 = plot(layout = (6, 5), size=(2400, 2000))
for (i, wait) in enumerate(0:5)
    if i <= length(empiricalMeans)
        difference = empiricalMeans[i]' - theoreticalMeans[i]'
        conditioned_difference = copy(difference)
        
        # Create a mask for the boolean condition, ensuring it matches the size of conditioned_difference
        mask = zeros(Bool, size(conditioned_difference))
        if eq3 == true
            mask[2:end, 1:end-1] = boolean_matrices[i] .== 1
            booltitle = "Eq3 <= 0"
        else
            mask = matricesPandemicEnded[i] .== 1
            booltitle = "Pandemic Ended"
        end
        # mask = matricesPandemicEnded[i] .== 1
        conditioned_difference[mask] .= 0

        minclim = min(minimum(empiricalMeans[i]), minimum(theoreticalMeans[i]))
        maxclim = max(maximum(empiricalMeans[i]), maximum(theoreticalMeans[i]))
        clims = (minclim, maxclim)
        
        # Only show x-labels for the last row
        xlabel_text = i == 6 ? "False Negative Rate" : ""
        
        heatmap!(p6[i, 1], false_negative_rates, trace_effeciency_rates, empiricalMeans[i]', 
                title="Empirical Mean (Wait=$wait)", xlabel=xlabel_text, 
                ylabel="Trace Efficiency", colorbar_title="Empirical Mean", clims=clims)
        
        heatmap!(p6[i, 2], false_negative_rates, trace_effeciency_rates, theoreticalMeans[i]', 
                title="Theoretical Mean (Wait=$wait)", xlabel=xlabel_text, 
                colorbar_title="Theoretical Mean", clims=clims)
        
        heatmap!(p6[i, 3], false_negative_rates, trace_effeciency_rates, 
                difference, 
                title="Difference Emp-Theo (Wait=$wait)", xlabel=xlabel_text, 
                 colorbar_title="Difference")
        heatmap!(p6[i, 4], false_negative_rates, trace_effeciency_rates, mask,
                title=booltitle, xlabel=xlabel_text, 
                colorbar_title="Pandemic Ended")
        heatmap!(p6[i, 5], false_negative_rates, trace_effeciency_rates, conditioned_difference,
                title="Difference for non-ended pandemics (Wait=$wait)", xlabel=xlabel_text, 
                colorbar_title="Empirical Trace Efficiency")
    end
end

display(p6)
savefig(p6, figurePath*"empirical_vs_theoretical_all_waits.png")



# With sampled Peff traced

empiricalMeans = []
theoreticalMeans = []

for wait in 0:5
    filename_emp = "code/OutputsMorten/empiricalMean_wait$(wait)NumExp20withSampledPefftraced.txt"
    filename_theo = "code/OutputsMorten/theoreticalMean_wait$(wait)NumExp20withSampledPefftraced.txt"
    
    if isfile(filename_emp) && isfile(filename_theo)
        empiricalMean = readdlm(filename_emp)
        theoreticalMean = readdlm(filename_theo)
        push!(empiricalMeans, empiricalMean)
        push!(theoreticalMeans, theoreticalMean)
        println("Loaded empirical and theoretical means for wait time $wait")
    else
        println("Files not found for wait time $wait")
    end
end

false_negative_rates = 0.0:0.1:1.0
trace_effeciency_rates = 0.0:0.1:1.0

# Create plot with 6 rows (one for each wait time) and 3 columns
using Plots


p7 = plot(layout = (6, 3), size=(1600, 1200))
for (i, wait) in enumerate(0:5)
        difference = empiricalMeans[i]' - theoreticalMeans[i]'
        minclim = min(minimum(empiricalMeans[i]), minimum(theoreticalMeans[i]))
        maxclim = max(maximum(empiricalMeans[i]), maximum(theoreticalMeans[i]))
        clims = (minclim, maxclim)
        
        # Only show x-labels for the last row
        xlabel_text = i == 6 ? "False Negative Rate" : ""
        
        heatmap!(p7[i, 1], false_negative_rates, trace_effeciency_rates, empiricalMeans[i]', 
                title="Empirical Mean (Wait=$wait)", xlabel=xlabel_text, 
                ylabel="Trace Efficiency", clims=clims)
        
        heatmap!(p7[i, 2], false_negative_rates, trace_effeciency_rates, theoreticalMeans[i]', 
                 ylabel="Trace Efficiency",title="Theoretical Mean (Wait=$wait)", xlabel=xlabel_text, 
                 clims=clims)
        
        heatmap!(p7[i, 3], false_negative_rates, trace_effeciency_rates, 
                difference, 
                 ylabel="Trace Efficiency",title="Difference Emp-Theo (Wait=$wait)", xlabel=xlabel_text 
                 )
end
display(p7)
savefig(p7, figurePath*"empirical_vs_theoretical_all_waits_sampledPtrace.png")