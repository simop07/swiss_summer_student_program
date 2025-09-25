#!/bin/bash
root -l -e '.L waveformAnalysisNeg.cpp+' \
        -e '.L waveformOsc.cpp+' \
        -e 'waveformAnalysis()' \
        -e 'waveformTotal()'