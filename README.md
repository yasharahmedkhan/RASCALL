### Installing the RASCALL package locally

1. Clone the repository and go into it. 

2. Activate the python environment by running 
```bash 
source pythonenv/bin/activate
```

3. Installing Dependencies  
You might run into errors or issues while installing dependencies. In that case, install each library needed manually with `pip install`. 
```sh 
pip install RASCALL 
``` 
and 
```sh
pip install -r requirements.txt
```

The tkinter library needs system level interface so install for linux with 
```sh
sudo apt-get install python3-tk
```
or follow the online doc on [installing tk](https://tkdocs.com/tutorial/install.html).

4. Some example commands
```sh
./rascall_plot --mol 'O=P(O)(O)O'
```

Other molecule samples can be found in `RASCALL_Molecule_Identifiers.xlsx`. Use their SMILE name on the left column. 






### Summarize Lorenzo Pico Ciamarra's work on issue finding code 

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