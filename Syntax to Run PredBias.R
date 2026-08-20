###This syntax file will run the PredBias Function.
###You first need to run the PredBias Function.R syntax
###to load the PredBias function into R's memory.


####Load data
##Set working directory (file location)
###Insert the directory you are using into the portion in quotes below.
###The input datafile should be in this directory and the results output file 
###will also be saved here

setwd("C:/PredBias")

###The input data file should be saved as a .csv
###comma-separated file.  Place the name of the
###input data file into the quotes below
###You can leave the R data file name as CombinedJeff
###or rename it to something else

PredBiasData<-read_csv("SampleData.csv")

###This statement runs the PredBias Function.
###You should modify the text as necessary
###The first portion in parentheses is the R input data
###file specified above.  The next two portions are the 
###variable names for the test and the criterion.  The last
###portion is the name of the output file.


PredBiasOutput<-PredBias(PredBiasData, "test", "criterion", "SampleOutput.csv", 
                         "SampleScatterPlot.png", "SamplePredictorHistogram.png",
                         "SampleCriterionHistogram.png")
