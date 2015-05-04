#!/usr/bin/env python
""" A small script to run the different algorithms of FitCrossCheckForLimits in parallel

Author: Nicolas Morange
Date:   2013-01-17
Email:  nicolas.morange@cern.ch

Description:
    This script is a simple wrapper to FitCrossCheckForLimits.
    The user sets the standard input and output variables, as well as the number of jobs
    to be launched in parallel.
    The user then indicates which algorithms to execute, with their parameters. For technical
    reasons some parameters to be provided here are dummy.
    Read FitCrossCheckForLimits documentation to find which ones are relevant to each algorithm.
    The script then executes the different algorithms in parallel, and merge the outputs

"""

import subprocess
import os
import shutil
import time

from ROOT import gROOT
from ROOT import TFileMerger



###############################################################
###################  User configuration  ######################
###############################################################

NCORES          = 5 # maximum number of jobs to launch in parallel
workspace       = "../../share/MCwithSignal/results/ws_LFV_combined_AllSYS_model.root" 
outputdir       = "results/MCwithSignal/"
workspaceName   = "combined"
modelConfigName = "ModelConfig"
ObsDataName     = "obsData"

# syntax: ["AlgName", mu, sigma, "IsConditional"]
# some parameters are not used by some algorithms. See the C++ clas for description.
# mu and sigma have different meanings for PlotsStatisticalTest
algs_to_execute = [
    # -----------------------------------------------------------------------------------
    # - Plot nominal and +/- Nsigma (for each nuisance paramater) for Data, signal+bkg
    # -----------------------------------------------------------------------------------
    ["PlotHistosBeforeFit",0.0,1.0,"false"],

    # -----------------------------------------------------------------------------------
    # - Control plots for morphing (ie, -1/0/+1 sigma --> continuous NP)
    # -----------------------------------------------------------------------------------
    ["PlotMorphingControlPlots",0.0,1.0,"false"],

    # ----------------------------------------------------------------------------------
    # - Plot histograms after unconditional fit (theta and mu fitted at the same time)
    # ----------------------------------------------------------------------------------
    ["PlotHistosAfterFitEachSubChannel",0.0, 1.0, "false"],
    ["PlotHistosAfterFitGlobal",0.0, 1.0, "false"],

    # -----------------------------------------------------------------------------------------
    # - Plot the conditionnal fitted nuisance parameters value (theta fitted while mu is fixed)
    # -----------------------------------------------------------------------------------------
    ["PlotHistosAfterFitEachSubChannel", 0.0, 1.0, "true"],
    ["PlotHistosAfterFitGlobal",0.0, 1.0, "true"],

    # -------------------------------------------
    # - Plot the nuisance parameters versus mu
    # -------------------------------------------
    ["PlotsNuisanceParametersVSmu",0.0, 1.0, "false"],

    # -------------------------------------------
    # - Plot the pulls and stat test from toys
    # -------------------------------------------
    ["PlotsStatisticalTest",1.0, 0.0, "false"],

]


###############################################################
################  End of User configuration  ##################
###############################################################


def main():
    """Parallelize the execution of algorithms in FitCrossCheckForLimits
    """

    print NCORES
    # cleaning
    outdir=outputdir.rstrip('/')
    outdir += '/'
    try:
        os.makedirs(outdir)
    except:
        pass

    # first, compile
    compile()

    # then, execute the different algorithms
    pids=[]
    logfiles=[]
    directories=[]

    for i,alg in enumerate(algs_to_execute):
        print i, alg
        if len(pids) >= NCORES: # manage number of jobs running
            wait_completion(pids)
        print "Launching job",i,":",alg
        output_f=open(outdir+"/output_"+str(i)+".log", 'w')
        logfiles.append(output_f)
        directory = outdir + str(i)
        directories.append(directory)
        formatted_args = alg[0] + "," + str(alg[1]) + "," + str(alg[2]) + "," + alg[3] \
                              + ",\"" + workspace + "\",\"" + directory + "\",\"" \
                              + workspaceName + "\",\"" + modelConfigName + "\",\"" + ObsDataName + "\""

        pids.append(subprocess.Popen(["root", "-l", "-b", "-q",
                                      "FitCrossCheckForLimits.C+("+formatted_args+")"],
                                     stderr=output_f, stdout = output_f))

    # Now just wait for completion of all jobs
    wait_all(pids)
    for f in logfiles:
        f.close()

    # and merge outputs
    print "Merging results..."
    tfm = TFileMerger()
    for directory in directories:
        tfm.AddFile(directory+"/FitCrossChecks.root")
    tfm.OutputFile(outdir+"/FitCrossChecks.root")
    tfm.Merge()

    for directory in directories:
        # TODO: pythonize the calls (shutil + os.path)
        subprocess.check_call(["cp","-r",directory+"/LatexFileNPs", outdir])
        subprocess.check_call(["cp","-r",directory+"/TextFileFitResult", outdir])
    print "Merging done !"

    # finally, remove the splitted files
    for directory in directories:
        shutil.rmtree(directory)
    print "All OK !"

def wait_all(pids):
    """Wait until completion of all launched jobs"""
    while len(pids)>0:
        wait_completion(pids)
    print "All jobs finished !"

def wait_completion(pids):
    """Wait until completion of one of the launched jobs"""
    while True:
        for pid in pids:
            if pid.poll() is not None:
                print "Process", pid.pid, "has completed"
                pids.remove(pid)
                return
        print "Waiting completion of jobs..."
        time.sleep(15) # wait 15 seconds before retrying


def compile():
    """Compile ROOT class"""
    gROOT.ProcessLine(".L FitCrossCheckForLimits.C+")

if __name__ == "__main__":
    main()
