package main

import (
	"fmt"
	"strconv"
	"math/rand"
    "math"	
	"time"	
	"os"
	"log"
)

func append_to_file (filename string, text []string) {
    // If the file doesn't exist, create it, or append to the file
    f, err := os.OpenFile(filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
    if err != nil {
        log.Fatal(err)
	}
	sep := " "
	for _, line := range text {	
		if _, err := f.WriteString(line+sep); err != nil {
			log.Fatal(err)
		}

	}
	if _, err := f.WriteString("\n"); err != nil {
		log.Fatal(err)
	}	
	if err := f.Close(); err != nil {
		log.Fatal(err)
	}	
}

func draw_children(R0 float64,days_preinfectious int,infectious_ends int,offspring_distribution string, child_location string, child_poisson []float64) []int {
	var children_array []int
	if (offspring_distribution == "poisson") {
		var number_of_children int = draw_poisson(R0)
		children_array := make_children_array(number_of_children,days_preinfectious,infectious_ends, child_location,child_poisson)
		return children_array

	} else if (offspring_distribution == "geometric") {
		var number_of_children int = draw_geometric(R0)
		children_array := make_children_array(number_of_children,days_preinfectious,infectious_ends, child_location, child_poisson)
		return children_array

	}

	return children_array
}

func get_child_location_array(infectious_ends int,child_location string,child_location_parameters [1]float64) []float64 {
	child_location_array := make([]float64,infectious_ends)
	if (child_location == "poisson") {
		mean := child_location_parameters[0]
		var denominator float64
		for day:= 0 ; day < infectious_ends ; day ++ {
			child_location_array[day] += math.Pow(mean,float64(day))/factorial(float64(day))*math.Exp(-mean)
			denominator += math.Pow(mean,float64(day))/factorial(float64(day))*math.Exp(-mean)
		}

		for day:= 0 ; day < infectious_ends ; day ++ {
			child_location_array[day] /= denominator
		}

	}
	return child_location_array
}

func draw_poisson(R0 float64) int {
	random_float := draw_random_float()
	var found_children int = 0
	var p_count float64 = 0
	var n_children float64 = -1
	for (found_children<1) {
		n_children += 1
		p_count += (math.Pow(R0,n_children)/factorial(n_children))*math.Exp(-R0)
		if (p_count >=random_float) {
			found_children += 3
		}
		
	}
	return int(n_children)

}

func draw_geometric(R0 float64) int {
	random_float := draw_random_float()
	var found_children int = 0
	var p_count float64 = 0
	var n_children float64 = -1
	for (found_children<1) {
		n_children += 1
		p_count += math.Pow((R0/(1+R0)),n_children)/(R0+1)
		if (p_count >=random_float) {
			found_children += 3
		}
		
	}
	return int(n_children)	
}

func draw_random_float() float64 {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	return r1.Float64()
}
func draw_random_integer(maximum int) int {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	return r1.Intn(maximum)
}
func make_children_array(n_children int,days_preinfectious int,infectious_ends int, child_location string, child_location_array []float64) []int {
	//var children_array [n_children]int
	children_array := make([]int,n_children)
	var integer int
	for child := 0 ; child < n_children ;child ++ {
		if (child_location == "uniform") {
			integer = draw_random_integer(infectious_ends-days_preinfectious)
		} else {
			integer = -1
			random_float := draw_random_float()
			found := 0
			var sum_poisson float64
			for (found == 0) {
				integer += 1
				sum_poisson += child_location_array[integer]
				if (sum_poisson >= random_float) {
					found = 1
				}
			}
		}
		children_array[child] = integer + days_preinfectious
	}
	
	return children_array
}

func factorial(k float64) float64{
	Res := 1
	if ( k>1) {
		for num := 2 ; num < int(k+1) ; num ++ {
			Res *= num
		}
	}
	return float64(Res)
}


func main() {


	// Definitions
	// -----------

	const time_max int = 50 // Maximum days run.
	const I_start int = 100 // Infectious at start
	const N_MaxPeople = 2000000 // How many people will we maximally get?

	const N_experiments int = 10

	// Epidemiological details
	const days_preinfectious int = 3 // Number of days before infectious
	const days_presymptomatic int = 4 // Number of days infectious but no symptoms
	const days_symptomatic int = 4 // Number of days infectious and showing symptoms
	const infectious_ends int = days_preinfectious+days_presymptomatic+days_symptomatic // First day after infectious (counting day 0 as preinfectious day)

	const days_between_tests int = 2 // This does not matter when people are never tested twice.


	const asymptomatic_fraction float64 = 0.3 // Fraction of infected that never get symptoms.

	const R0 float64 = 2.0//3//2.5 // Mean number of children in full period of infection.
	const offspring_distribution string = "poisson"
	//const offspring_distribution string = "geometric"

	
	//const child_location string = "uniform"
	const child_location string = "poisson"
	child_location_parameters := [1]float64{3}
	child_location_array := get_child_location_array(infectious_ends,child_location,child_location_parameters)

	// Societal details
	var wait_test_taken int = 0  // Number of days before test is taken
	var wait_test_result int = 0 // Number of days before test result arrives after test is taken

	// Test-and-trace details
	//var false_negative_test_rate float64 = 0.00 // Probability that a test of a COVID-positive comes back negative
	var tracing_efficiency float64 = -0.02//+34*0.02 // Fraction of children that are found through contact tracing

	

	// --------------------




	
	var filename string = "TestSensitivity_Istart"+strconv.Itoa(I_start)+"Nexp"+strconv.Itoa(N_experiments)+"R0"+fmt.Sprintf("%f", R0)+"_Dayspreinfectious"+strconv.Itoa(days_preinfectious)+"Dayspreinfectious"+strconv.Itoa(days_presymptomatic)+"Dayssymptomatic"+strconv.Itoa(days_symptomatic)+"Testwait"+strconv.Itoa(wait_test_taken)+"Resultwait"+strconv.Itoa(wait_test_result)+"Asymptomatics"+fmt.Sprintf("%f", asymptomatic_fraction)+"ChildLocation"+child_location+"OffspringDist"+offspring_distribution+".txt"
	valuesText1 := []string{fmt.Sprintf("False negative test rate,"+"Tracing efficiency,"+"N_infected_done,"+"N_recovered,"+"Reff")}

	append_to_file("Outputs/"+filename,valuesText1)


	for d_trace :=0 ; d_trace < 51; d_trace ++ {
		tracing_efficiency += 0.02
		var false_negative_test_rate float64  = -0.02
		fmt.Println("Parameters:", "False neg:",false_negative_test_rate, "Trace efficiency:",tracing_efficiency)
		for d_false := 0 ; d_false < 51 ; d_false ++ {
			false_negative_test_rate += 0.02
			var R0_mean float64
			var N_recovered_mean int = 0 // How many have now recovered?
			var N_infected_done_mean int = 0 // How many have now recovered?			
			for exp:= 0 ; exp < N_experiments ; exp ++ {

				// Definitions: Keeping count of results
				var N_infected int = I_start+0 // Count how many have been infected
				var N_recovered int = 0 // How many have now recovered?
				var N_infected_done int = 0 // How many have now recovered?
				
				// Keeping count of states

				States := make([]int,N_MaxPeople) // Array with Day on entry. Negative if not infected,
				Test_waiting := make([]int,N_MaxPeople) // Array with Day on entry. Negative if not waiting,
				Result_waiting := make([]int,N_MaxPeople) // Array with Day on entry. Negative if not waiting,
				Start_tracing := make([]int,N_MaxPeople) // Array with Day on entry. Negative if not waiting,
				Test_twice := make([]int,N_MaxPeople) 
				
				Asymptomatic:= make([]int,N_MaxPeople) // Array with 0 or 1 on entry. 0 if normal, 1 if always asymptomatic.

				// Define states, waiting times of non-seeds. To be changed when infected, waiting for test, etc..
				for entry:=0 ; entry < N_MaxPeople ; entry ++ {
					if (entry >= N_infected) {
						States[entry] -= 9
					} 
					random_float := draw_random_float()
					
					// See if asymptomatic individual
					if (random_float<=asymptomatic_fraction) {
						Asymptomatic[entry] = 1 +0
					}
					
					Test_waiting[entry] -= 9
					Result_waiting[entry] -= 9

					

				}
				// Decide when seeds will infect other people. Saved as array with entries equal to the day on which a child is infected.
				Infect_children := make([][]int,N_MaxPeople)
				for entry:=0 ; entry < N_infected ; entry ++ {
					Infect_children[entry] = draw_children(R0,days_preinfectious,infectious_ends,offspring_distribution,child_location, child_location_array)
				}
				// Make array to keep track of a node's children. To be used in Contact tracing and when calculating Reff at last.
				Children := make([][]int,N_MaxPeople)


				for t := 0 ; t < time_max ; t++ {



					// Advance all infected 1 time step
					for entry:=0 ; entry < N_MaxPeople ; entry ++ {
						if (States[entry]>=0 && States[entry]!=infectious_ends) {
							States[entry]+= 1

							// If node ends its infectious period....
							if (States[entry]==infectious_ends) {
								N_recovered += 1 // One more infectious has infected everybody it will get to infect..
								N_infected_done += len(Children[entry]) // Add to the number of children 							
							}

						}
						// If infectious node is waiting for a test...
						if (Test_waiting[entry] >= 0) {
							// If test is taken today, start waiting for result
							if (Test_waiting[entry] == 0) {
								Result_waiting[entry] = wait_test_result//0 +0
								if (Test_twice[entry] == 1) {
									Test_waiting[entry] = days_between_tests +0
									Test_twice[entry] = 0
								}


							}
							Test_waiting[entry] -= 1 // advance wait 1 day.


						}
						// If infectious node is waiting for test result...
						if (Result_waiting[entry] >= 0) {
							// If test result comes back today..
							if (Result_waiting[entry] == 0) {

								// Draw float to determine whether test is false negative.
								random_float := draw_random_float()
								// Important: With the following definition, we can only test positive if people are infectious...
								if (States[entry] >= 0 && States[entry]-wait_test_result >= days_preinfectious && States[entry]-wait_test_result  < infectious_ends &&  random_float > false_negative_test_rate) {
									// Test comes back positive. 
									Start_tracing[entry] = 1 +0 // Start tracing
									States[entry] = infectious_ends // Place node in quarantine

									// Now cancel future tests / results.
									Test_twice[entry] = 0
									Test_waiting[entry] = -9
									Result_waiting[entry] = -9
									
									// Now keep count of how many the node infected during its infectious period.
									N_recovered += 1 // One more infectious has infected everybody it will get to infect..
									N_infected_done += len(Children[entry]) // Add to the number of children 
								} //else if (States[entry] >= 0 && States[entry]-wait_test_result-wait_test_taken >= days_preinfectious && States[entry]-wait_test_result-wait_test_taken  < infectious_ends &&  random_float <= false_negative_test_rate) {
									// Test comes back false negative. Reset the node's waiting time for test and result (in case nodes will try to get tested again later)
									//Test_waiting[entry] = -9
									//Result_waiting[entry] = -9
								//}
							}
							// Advance waiting if node is still waiting (could have stopped waiting if result was false negative)
							if (Result_waiting[entry]>=0) {
								Result_waiting[entry] -= 1
							}
						}
					}		

					// Infect all children on this timestep
					for entry:=0 ; entry < N_MaxPeople ; entry ++ {
						// If a node is infectious and not in quarantine, waiting for test
						if (States[entry]>=days_preinfectious && States[entry] < infectious_ends && Test_waiting[entry]<0 && Result_waiting[entry] < 0  ) {
							// Check if the node is supposed to get a child this time step
							for child:=0 ; child < len(Infect_children[entry]) ; child ++ {
								if (Infect_children[entry][child]==States[entry]) {
									// If node gets child, update state of child and when children will get infected
									States[N_infected] = 0 +0
									Infect_children[N_infected] = draw_children(R0,days_preinfectious,infectious_ends,offspring_distribution,child_location, child_location_array)
									Children[entry] = append(Children[entry],N_infected) // save child as child of node "entry"

									N_infected += 1


								}
							}
						}
					}

					// Get test if getting symptoms
					for entry:=0 ; entry < N_MaxPeople ; entry ++ {
						if (States[entry] == days_preinfectious+days_presymptomatic && Test_waiting[entry]<0 && Result_waiting[entry]<0 && Asymptomatic[entry] == 0) {
							Test_waiting[entry] = wait_test_taken +0
						}

						// Trace
						if (Start_tracing[entry] == 1) {
							for child := 0 ; child < len(Children[entry]) ; child ++ {
								random_float := draw_random_float()
								// Toss coin if child is going to be found.
								if (random_float<=tracing_efficiency) {
									// If successful, Test child.
									Test_waiting[Children[entry][child]] = int(math.Max(float64(3-States[Children[entry][child]]),float64(wait_test_taken)))//wait_test_taken +0
									//Test_twice[Children[entry][child]] = 1 // Test twice if found through tracing..
								}
							}

							Start_tracing[entry] = 0
						}

					}


					

				}
				R0_mean += (float64(N_infected_done)/float64(N_recovered))/float64(N_experiments)
				N_recovered_mean += N_recovered
				N_infected_done_mean += N_infected_done
			}


			valuesText := []string{fmt.Sprintf("%f", false_negative_test_rate)+","+fmt.Sprintf("%f", tracing_efficiency)+","+strconv.Itoa(int(N_infected_done_mean))+","+strconv.Itoa(int(N_recovered_mean))+","+fmt.Sprintf("%f", R0_mean)}

			append_to_file("Outputs/"+filename,valuesText)
			
		}
	}
}
