# Protocol_EMG

## REQUIREMENTS

MATLAB 2018 or higher versions

## INSTALL

Download the program and relative directories

# USAGE

Launch [Indexes, GaitPhases] = SbS_indexesFromEMGinEBformat(filenameIN, CHforIdxFilename, OUTPUTdir) 
[Output terms are optional]

Results will be saved authomatically in the output directory.

AIM

The program will load the FILENAME csv file with the EMG signals and calculate the activation indexes for each channels in the CHforINDEXES file.
The cocontraction index will be calculated for each couple in the CHforINDEXES file.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                                                                       %%%%%%%%%%%%
%%%%%%%%%%%%  [Indexes, GaitPhases] = SbS_indexesFromEMGinEBformat(filenameIN, ...                 %%%%%%%%%%%%
%%%%%%%%%%%%                          PHASESfilenameIN, CHforIdxFilename, OUTPUTdir)               %%%%%%%%%%%%
%%%%%%%%%%%%                                                                                       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                                                         %%%
%%%         Eurobench STEPbySTEP routine to calculate SbS indexes from EMG in Eurobench format              %%%
%%%                                                                                                         %%%
%%%         IN:     filenameIN          --->    STRING                                                      %%%
%%%                                             '.\INPUTdir\filename.csv' with EMG data in EB format        %%%
%%%                                                                                                         %%%
%%%                 PHASESfilenameIN   --->    STRING                                                       %%%
%%%                                             '.\INPUTdir\filename.ylm' with gait phases in EB format     %%%
%%%                                                                                                         %%%
%%%                 CHforIdxFilename    --->    STRING                                                      %%%
%%%                                             '.\INPUTdir\filename.yml' file with couples of numbers      %%%
%%%                                             to calculate cocontraction                                  %%%
%%%                                             e.g. CHforINDEXES: [ [ [6 7], [10 11], [15 16] ]            %%%
%%%													    %%%
%%%                 OUTPUTdir           --->    STRING                                                      %%%
%%%                                             '.\OUTPUTdir output directory where indexes are saved       %%% 
%%%                                                                                                         %%%
%%%                                                                                                         %%%
%%%         OUT:    Indexes             --->    STRUCT                                                      %%%
%%%                                             Fields:                                                     %%%
%%%                                                     ActivationIndexes(i)                                %%%
%%%                                                     	activation index of EMG i-channel               %%%
%%%                                                     CocontractionIndexes(j)                             %%%
%%%                                                         CocontractionIndexes EMG j-couple (e.g. [1 2])  %%%
%%%                 GaitPhases          --->    STRUCT                                                      %%%
%%%                                             Fields:                                                     %%%
%%%                                                     Time                                                %%%
%%%                                                     Perc                                                %%%
%%%                                                                                                         %%%
%%%         Example:    [Indexes, GaitPhases] =                                                             %%%
%%%                     SbS_EMG_indexes('.\input\input_emg_file.csv', '.\input\input_emg_file.csv', ...     %%%
%%%                     'input\CHforINDEXES.yml','.\output\'                                                %%%
%%%                                                                                                         %%%
%%%                                                                                                         %%%
%%%         SUBROUTINES:            ReadYaml, FiltButterLBH, EnvelopeHilbert, ...               	        %%%
%%%                                 discrete_integrate, cocontraction_winter, 				                %%%
%%%                                                                                                         %%%
%%%                                                                                                         %%%
%%%         Author:     Marco Caimmi                                                                        %%%
%%%                     STIIMA Nationl Research Council of Italy                                            %%%
%%%                     marco.caimmi@stiima.cnr.it                                                          %%%
%%%                     marco.caimmi@gmail.com                                                              %%%
%%%                                                                                                         %%%
%%%         Year:       2021                                                                                %%%
%%%                                                                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```


## ACKNOWLEDGEMENTS

<a href="http://eurobench2020.eu">
  <img src="http://eurobench2020.eu/wp-content/uploads/2018/06/cropped-logoweb.png"
       alt="rosin_logo" height="60" >
</a>

Supported by Eurobench - the European robotic platform for bipedal locomotion benchmarking.
More information: [Eurobench website][eurobench_website]

<img src="http://eurobench2020.eu/wp-content/uploads/2018/02/euflag.png"
     alt="eu_flag" width="100" align="left" >

This project has received funding from the European Union’s Horizon 2020
research and innovation programme under grant agreement no. 779963.

The opinions and arguments expressed reflect only the author‘s view and
reflect in no way the European Commission‘s opinions.
The European Commission is not responsible for any use that may be made
of the information it contains.

[eurobench_logo]: http://eurobench2020.eu/wp-content/uploads/2018/06/cropped-logoweb.png
[eurobench_website]: http://eurobench2020.eu "Go to website"
