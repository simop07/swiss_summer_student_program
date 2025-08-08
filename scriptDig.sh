#!/bin/bash
root -l -e '.L waveformAnalysisPos.cpp+' \
        -e '.L waveform.cpp+' \
        -e 'waveformAnalysis()' \
        -e 'waveformTotal()'