#!/bin/bash
root -l -e '.L waveformAnalysisNeg.cpp+' \
        -e '.L waveformOscDownsample.cpp+' \
        -e 'waveformAnalysis()'