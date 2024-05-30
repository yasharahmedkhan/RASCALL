
# RASCALL 1.0
The main function of RASCALL is that given inputs of a molecular dictionary and a functional group dictionary, output approximate spectral data for molecules in the gas phase. For more see the [RASCALL paper](https://pubs.rsc.org/en/content/articlehtml/2019/cp/c8cp07057a). 


# Quick Start

### Installing the RASCALL package locally
1. Clone the repository and go into it. 

2. Activate the python environment by running 
```bash 
source pythonenv/bin/activate
```

3. Installing Dependencies
   
In order to install dependencies, once you are in the virtual environment, run the shell script titled "setup_rascall.sh".
To do this, run the following commands:
```sh 
chmod +x setup_rascall.sh 
``` 
then 
```sh
./setup_rascall.sh
```
This should check and install all the dependencies required in the virtual environment, and allow you to use RASCALL effectively. This is in place of the previous "requirements.txt" installation that was causing some issues.  

4. Some example commands to test the installation
```sh
./rascall_plot --mol 'O=P(O)(O)O'
```
Plots a functional group, one plot per bond in the database. Plots should appear in a separate window.  

```sh
./rascall_list --mol 'O=P(O)(O)O'
```
Plots a functional group, one plot per bond in the database. Plots should appear in a separate window.  

There are more commandline options, see the code files to examine all the command line options available. Currently, the useful ones are: 
- `--fg=[functional group]` for instance `--fg=[!#1]C#C[!#1]` 
- `--mol=[molecoleSMILEname]` for example  `--mol='CN(CC)C'` 
- `--mf=["halo","hydro","all"]` 

Other molecule samples and their functional groups can be found in `RASCALL_Molecule_Identifiers.xlsx`. Use their SMILE name on the left column. 



# Dev & Contributing 

## Getting started

Skim the [RASCALL paper](https://pubs.rsc.org/en/content/articlehtml/2019/cp/c8cp07057a) which describes the function, workings and science of RASCALL. Section 5.3 `Future versions` includes plans for updates on code base. The RASCALL paper and this repository will serve as the main reference. Notes by previous contributors can be found in the `notes/` folder. 

The quickest way to understand the workings of RASCALL is to run the code in Quick Start, read the code, make small modifications and see what they are doing. 


## Intermediate Projects
- fix `rascall_plot` line 21 comment on frequency window



## Lorenzo Pico Ciamarra's work on issue finding code
Lorenzo's detailed notes on issue_finding is in ```notes/LorenzoRASCALL_report.docx```.

To use his code, run  
```bash 
python3 executable_issues_finder.py
```

From Lorenzo: 
I am also including here a zipped folder with sample text file outputs for a few functional groups, which needs to be unzipped (but kept within a "func_group_stats" folder) and placed directly inside the RASCALL master directory (without changing any names). With this, using the plot-only function I mention in my report, you should be able to reproduce also the histogram plots which form the final output of the code. Of course, I'm happy to answer any questions or doubts you might have on this code!

autofit.py
funcs_histograms.py
executable_issues_finder.py
NIST_funcs.py
functional_analysis.py
plot_histogram.py
functional_parser.py
func_group_stats.zip
RASCALL_report.docx

I just noticed that, in fact, I did not explain in my manual how to use this to only plot the histograms without performing new calculations: after having placed all the files in the correct directories as explained in the manual, just run exectuable_issues_finder.py and, when prompted to do so at the very beginning, type "s" in the terminal (which stands for skipping the calculations and going directly to the plotting). It will then ask which kinds of histograms you would like to plot - I explain the significance of each of them in the report. Hope that's clear!


### Improvements on Lorenzo's work 
- ```autofit.py``` that estimates the spectrum can be improved in terms of efficiency and accuracy. 
- Currently the path for `executable_issues_finder.py` is designed for Windows operating system. The code will need to be modified to be used for other OS.   
