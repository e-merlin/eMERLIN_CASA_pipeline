
## TL;DR ##
This is the CASA pipeline for eMERLIN data, it is designed to be fully parallelised. It is, however, very early stages and will progress rapidly. It currently does work but does not do calibration yet. Feel free to change/modify.

## Dependencies ##
- aoflagger v2.9 (see https://sourceforge.net/projects/aoflagger/)
- CASA v4.6+ (see https://casa.nrao.edu/), tested on CASA v4.6.0
- python2.7
- git

## Usage ##
The pipeline uses a GUI interface as an input and now has a headless version. To run, do the following:

1. Make a new directory
2. cd into directory and mkdir data
3. `git clone https://github.com/e-merlin/CASA_eMERLIN_pipeline.git`
4. Copy data (either fits or ms) into data
5. Run the pipeline
  * If using the input file/ headless version (input.txt included in CASA_eMERLIN_pipeline) run:
    a. cp inputs.txt to cwd and edit as appropriate
    b. `mpicasa -n num_cores casa -c CASA_eMERLIN_pipeline/eMERLIN_CASA_pipeline.py -i inputs.txt`
  * Else to run GUI
    `mpicasa -n num_cores casa -c CASA_eMERLIN_pipeline/eMERLIN_CASA_pipeline.py -g`

To run not in parallel remove `mpicasa -n num_cores`

IMPORTANT: The pipeline will look in your data directory for a ms or mms with the prefix defined in inputs and copy it to your current directory if it cannot find it in the current directory.

N.b. replace num_cores with the number of cores you would like to use
