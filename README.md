# Protocol_EMG

## Aim
Code to calculate the EMG performance indicators (activation and cocontraction indexes) whithin the STEPbySTEP project

## Requirements

MATLAB 2018 or higher versions

## Install

Download the program and relative directories
Insert the YAML directory in the Matlab path

## Functioning

Launch [Indexes, GaitPhases] = SbS_EMG_indexes(FILENAME, CHforIdxFilename, Nexus, toBePlot)

(The output terms are optional)

Results will be saved authomatically in the output directory.


The program will load the FILENAME csv file with the EMG signals and calculate the activation indexes for each channels in the CHforINDEXES file.
The cocontraction index will be calculated for each couple in the CHforINDEXES file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%                                                                                       %%%%%%%%%%%%

%%%%%%%%%%%% [Indexes, GaitPhases] = SbS_EMG_indexes(FILENAME, CHforIdxFilename, Nexus, toBePlot)  %%%%%%%%%%%%

%%%%%%%%%%%%                                                                                       %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%                                                                                                         %%%

%%%         Eurobench STEPbySTEP routine to calculate SbS EMG indexes                                       %%%

%%%                                                                                                         %%%

%%%         IN:     filename            --->    STRING                                                      %%%

%%%                                             name of the CSV file with EMG and                           %%%

%%%                                             acceleration data from the acquisition system               %%%

%%%                                                                                                         %%%

%%%                 CHforIdxFilename    --->    YAML file with couples of numbers                           %%%

%%%                                             CHforINDEXES to calculate cocontraction                     %%%

%%%                                             e.g. CHforINDEXES:                                          %%%

%%%                                                                 - [1 2]                                 %%%

%%%                                                                 - [1 3]                                 %%%

%%%                                                                 - [2 3]                                 %%%

%%%                                                                                                         %%%

%%%                 Nexus               --->    BOOL                                                        %%%

%%%                                             0 acquired with Delsys software                             %%%

%%%                                             1 acquired with Nexus software                              %%%

%%%                                                                                                         %%%

%%%                 toBePlot            --->    BOOL                                                        %%%

%%%                                             0 do not plot the data                                      %%%


%%%                                             1 plot the data                                             %%%

%%%                                                                                                         %%%

%%%         OUT:    Indexes             --->    STRUCT                                                      %%%

%%%                                             Fields:                                                     %%%

%%%                                                     ActivationIndexes                                   %%%

%%%                                                     CocontractionIndexes                                %%%

%%%                 GaitPhases          --->    STRUCT                                                      %%%

%%%                                             Fields:                                                     %%%

%%%                                                     Time                                                %%%

%%%                                                     Perc                                                %%%

%%%                                                                                                         %%%

%%%                                                                                                         %%%

%%%                                                                                                         %%%

%%%                 SUBROUTINES:    ReadYaml,LoadDelsysData, FiltButterLBH, EnvelopeHilbert, ...            %%%

%%%                                 discrete_integrate, cocontraction_winter, WriteYaml                     %%%

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


