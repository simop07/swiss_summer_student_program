#!/bin/bash
root -l -e '.L waveformAnalysisNeg.cpp+' \
        -e '.L waveformOscDownsample30.cpp+' \
        -e 'waveformAnalysis()' \
        -e 'waveformTotal()'