# AnalysisVisualization
#Introduction and Concept
Back-end for visualizing data in analyses using ROOT files with an arbitrary cut-string interpreter.
In principle the skeleton of what does selection on a dataframe of any kind is comprised of the classes within the utilities folder
"utils/Parametrization.py". Parametrization.py contains the concept of a category and process.
Categories are selections on the dataset of interest with certain cuts ... and may contain sub-categories as well.
Processes make up the raw distributions that will be considered - for most cases this is data and Monte Carlo (MC).
To utilize these classes, a main script is made to interface with them on the analysis of interest. This main script contains the "cut string interpreter".
yis the main script that does all the heavy lifting and contains the cut-string interpreter. The classes are invoked within the function and can be interfaced to the analyzer's content with yaml files.
The recommended way to use these classes and this script is to make a *yaml* file containing the categories and another separate *yaml* file for all the distributions or processes the analyzer should consider.
Conveniently, with the ZH_Run2 framework - a csv file containing all the MC distributions is listed. Therefore the current implementation just reads this in and creates all the processes of interest except for those pertaining to data (data itself and data-driven methods)
In addition, it may be useful to create new variables on the fly or cut on functions of a given variable. These are both done in the framework.
To create a new variable, add it to the new variables section within the category example - this interfaces with the *functions.py* script in the utils folder to create a new numpy array from existing variables.
If implemented in the same way as the **MakeDistributions_HAA_2016.py** script the new variables can even be used in cuts on the fly.

# Setup
Setup the right scram architecture and clone then source the setup:
```
scram pro -n "nanov7_basic_10_6_4" CMSSW_10_6_4
cd nanov7_basic_10_6_4/src  
cmsenv
git clone https://github.com/samhiggie/AnalysisVisualization.git
cd AnalysisVisualization
bash setup.sh

```

# Plotting
The main script that controls plotting from skimmed NanoAOD is **MakeDistributions_HAA_2016.py**
To make ultra skimmed ntuples this script can be used with just MC by invoking (*-s* for skimmed vs histogram generation)
```
python MakeDistributions_HAA_2016.py -c cat_mmtt_2016.yaml  -csv MCsamples_2016_v6_yaml.csv  -i /eos/home-s/HAA_ntuples/June2020/ -p processes_special_mmtt.yaml -o main_outputdirectory -fo fake_factor_measure_outputdirectory -ch mmtt -s
```
To invoke only data-driven measurement (notice the *dmZH* flag) done with histograms!
```
python MakeDistributions_HAA_2016.py -c cat_mmtt_2016.yaml  -csv MCsamples_2016_v6_yaml.csv  -i /eos/home-s/HAA_ntuples/June2020/ -p processes_special_mmtt.yaml -dmZH-o main_outputdirectory -fo fake_factor_measure_outputdirectory -ch mmtt
```
To conduct the data-driven application (once the measurement is done and the fake factor histograms exist) - note this is to produce ultra skimmed
```
python MakeDistributions_HAA_2016.py -c cat_mmtt_2016.yaml -csv MCsamples_2016_v6_yaml.csv -i /eos/home-s/HAA_ntuples/June2020/ -p processes_special_mmtt.yaml -o main_output_file -ch mmtt -s -ddZH -fi fake_factor_measure_outputdirectory
```
After the ultra skimmed files are made one can then plot the histograms of associated variables in a short amount of time via (for datadriven use *-ddZH* flag)
```
python MakePlots_skimmed.py -i skimmed_mainoutputfile_mmtt_inclusive.root -o main_outputdirectory -c cat_mmtt_2016.yaml --ch mmtt -ddZH
```

# Automation
A bash script exists to conduct all the previous commands in an automated way for a single run per channel.
```
nohup bash runFF_mmtt.sh > mmtt.out &
```

# Output
Several features are implemented for checking the output. The main thing
