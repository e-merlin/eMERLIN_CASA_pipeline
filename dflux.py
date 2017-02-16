import numpy
import math
import sys


#######
# Python version of 3C286 flux calculation program (original author unknown)
# ..............................
# Author DMF       20/10/2011
# ..............................
#
# Update to use Perley & Butler 2012 coefficients
# 10/04/2013
# DMF
########

# Reworked to use the 1999 VLA flux formula, and a 2nd formula to give a continuous estimate of the resolved fraction, by Ian Stewart, JBO, 8 Aug 2007.
# Minor changes by amsr, 8 Aug 2007

# my $program_name = 'dflux'; # $0 returns the './' prefix if this is used.


def dfluxpy(freq,baseline):


    lowest_freq = 300.0;
    highest_freq = 50000.0;
    if (freq < lowest_freq or freq > highest_freq):
        print "Frequency must be between $lowest_freq and $highest_freq MHz. \n"

    # Old values for 3C286
    # A = 1.23734
    # B = -0.43276
    # C = -0.14223
    # D = 0.00345

    # New values taken from AIPS SETJY 31DEC11
    # Values as of 2010

    # A = 1.2361
    # B = -0.4127
    # C = -0.1864
    # D = 0.0294
    
    # Perley & Butler 2012 values
    A = 1.2515
    B = -0.4605
    C = -0.1715
    D = 0.0336

    log10f = (math.log(freq)/2.3025851) - 3.0; # Why the -3? Because freq has to be GHz for the formula to work.
    log_flux = A + B*log10f + C*log10f*log10f + D*log10f*log10f*log10f
    vlaflux = math.pow(10.0,log_flux)


    

    # The VLA flux must now be corrected to account for the higher resolving power of merlin. The formula used was obtained with the help of Peter Thomasson. If we assume that 3C286 is represented by a gaussian of angular size theta_s, and represent the resolving power of the given baseline as a function of frequency f and antenna separation L by theta_b(f,L), then the reduction in central flux density A(0) due to the finite theta_s is given by
    #
    #	                           1
    #	            -----------------------------------
    #	 A'(0)       2 pi (theta_b(f,L)^2 + theta_s^2)
    #	------- = --------------------------------------- ,
    #	 A(0)                      1
    #	                   ---------------------
    #	                    2 pi theta_b(f,L)^2
    #
    #	               1
    #	        = -------------- ,
    #	           1 + rho(f,L)
    #
    # where the resolved fraction rho(f,L) is given by
    #
    #	              theta_s^2
    #	rho(f,L) = ---------------- .
    #	            theta_b(f,L)^2
    #
    # Use of theta_b(f,L) = k/(fL) allows this to be written
    #
    #	           (   f*L     )^2
    #	rho(f,L) = (-----------)   * rho_ref .
    #	           (f_ref*L_ref)
    #
    # The reference value of rho is fixed at 0.04 for the MK-TA baseline at 5 GHz (Peter Thomasson).

    ref_bl_length = 11236.79 # MK-TA separation in metres.
    ref_freq = 5000.0
    ref_rho = 0.04
    thisbl = "this baseline (Mk-Ta)"
    
    bl_length = baseline
    # bl_str = sprintf "%8.2f", $ref_bl_length;

    
    frac = (freq / ref_freq) * (bl_length / ref_bl_length)
    rho = frac * frac * ref_rho
    merlinflux = vlaflux / (1.0 + rho)

    # Another useful quantity is the resolved percentage:
    #
    resolved_percent = 100.0 * rho / (1.0 + rho)
    caution_res_pc = 10.0

    return vlaflux, merlinflux, resolved_percent, caution_res_pc, thisbl





#####################################################################
'''

if __name__ == '__main__':

    helptext = "This command prints the flux density of the flux calibrator 3C286\nat a given frequency along with advice for use with MERLIN. \n\nUsage: dfluxpy [-h]  <freq (MHz)> [<baseline (m)>]\n Examples:\n dflupy 4994 assumes MK-TA baseline.\n dfluxpy 6034.8 15923 ie, work out the MERLIN correction for the DA-TA baseline instead.\n Some short MERLIN baselines (m):\n  MK-TA: 11236.79 \n  DA-TA: 15923.55\n  MK-DA: 17737.45\n  CA-MK is about 200000 m depending on hour angle (not recommended).\n"


    if (len(sys.argv)<2):
        print helptext
        # raise 'notEnoughInfoOnCmdLine'
        sys.exit()



    elif (len(sys.argv)>2):
        bl_length = float(sys.argv[2])
        if (bl_length <= 0.0):
            print "Second command-line argument does not seem to be non-zero baseline distance."
        thisbl = "this baseline"
        freq = float(sys.argv[1])

    else:
        freq = float(sys.argv[1])
        bl_length = 11236.79


    print freq, bl_length
    vlaflux, merlinflux, resolved_percent, caution_res_pc, thisbl  = dfluxpy(freq,bl_length)
    


    print '\nAt frequency',freq, 'MHz:'
    print 'Flux of 3C286 =',vlaflux, 'Jy using the VLA formula with 1999.2 coefficients.'
    print '3C286 is resolved by',resolved_percent, 'percent on a baseline length of ',str(bl_length),' m.'
    if (resolved_percent >= caution_res_pc):
        print' ->(CAUTION! Resolved values higher than about',caution_res_pc, 'percent become increasingly inaccurate.)'
    print 'Suggested flux density for MERLIN on', thisbl, 'is' ,merlinflux, 'Jy. \n'

'''
