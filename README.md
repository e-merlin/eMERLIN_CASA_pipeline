#CASA eMERLIN Pipeline

## TL;DR ##
This is the CASA pipeline for eMERLIN data, it is designed to be fully parallelised. It is, however, very early stages and will progress rapidly. It currently does work but does not do calibration yet. Feel free to change/modify. 

## Dependencies ##
- aoflagger v2.9 (see https://sourceforge.net/projects/aoflagger/)
- CASA v4.7.1 (see https://casa.nrao.edu/)
- python2.7 

## Usage ##
The pipeline only works in parallel mode at the moment and uses a GUI interface as an input. To run, move to your current working directory with the UVFITS (.fits), measurement set (.ms) or multi-measurement set (.mms) file, extract the pipeline to that directory and run: 

`mpicasa -n num_cores casa -c CASA_eMERLIN_pipeline.py`

N.b. replace num_cores with the number of cores you would like to use



