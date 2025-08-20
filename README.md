# UZH Summer Student Program

[![DOI](https://zenodo.org/badge/1024295677.svg)](https://doi.org/10.5281/zenodo.16323819)
 [![GitHub 
license](https://img.shields.io/github/license/simop07/swiss_summer_student_program)](https://github.com/simop07/swiss_summer_student_program/blob/main/LICENSE)

## Repository contents (work in progress)
This repository contains the material developed during my time as Visiting Student at the [University of Zurich](https://www.uzh.ch/en.html) (UZH) as part of the [Swiss Summer Student Program](https://swiss.sspp.program.phys.ethz.ch/) in Particle Physics, held from 14 July to 14 September 2025.

## Introduction
The aim of my project is to analyse the reflectance and absorption of light in Polytetrafluoroethylene (PTFE), often referred to as _Teflon_ by the brand name. This material is commonly used in the construction of liquid xenon detectors due to its high diffuse reflectivity for VUV (Vacuum Ultra-Violet) scintillation light, a type of radiation primarily produced as a consequence of particle's interaction with noble gases. As liquid xenon detectors have proved effective in the search for rare processes, such as elastic scattering of dark matter particles off nuclei, they are largely utilised in dark matter research for their great sensitiviness and performance.

To efficiently collect the scintillation light produced by particles interacting with liquid xenon, VUV reflectors are required. While PTFE is excellent for this purpose, it is desired to decrease its overall amount in order to minimise the radiogenic background produced by its inherent impurity, without loosing light collection efficiency. Hence, the complete knowledge of the scintillation light transmitted in PTFE reveals fundamental, expecially in low-background experiments.

The first part of my project is devoted to measure the reflectance and absorption of light in PTFEs by measuring the scintillation light detected by PMTs when light is emitted from a blue LED. To analyse the signals measured by the PMTs, a waveform and pulse analysis is performed through a software based on `ROOT` framework and `C++`.

### Pulse analysis
- [merge.cpp](merge.cpp) is just used to merge multiple waveforms (one per .scv or .xlsx file) into one single .txt file.

- [compare.cpp](compare.cpp) is used to compare 1D histos retrieved from the analysis.

- [lightAnalysis.cpp](lightAnalysis.cpp) is used to compute transmittance and reflectance in various time regions.

To perform pulse analysis, two lines of development are followed depending on the used logic:
- [waveformAnalysisPos.cpp](waveformAnalysisPos.cpp), pulses with positive logic;
- [waveformAnalysisNeg.cpp](waveformAnalysisNeg.cpp), for pulses with negative logic.

In addition, two main source files are presented depending on the used data:
- [waveform.cpp](waveform.cpp), for digitised data;
- [waveformOsc.cpp](waveformOsc.cpp), for data from the oscilloscope.
- [waveformOscDownsample.cpp](waveformOscDownsample.cpp), for downsampled data from the oscilloscope.

The aim of pulse analysis is to extrapolate several parameters of interest which will leter be used to find the reflectance and transmittance of light in PTFEs. You can see the definition of parameters of interest directly in [waveform.cpp](waveform.cpp) and [waveformOsc.cpp](waveformOsc.cpp).

Graphs in [pulses.pdf](plots/pulses.pdf) and [pulses_osc.pdf](plots/pulses_osc.pdf) show pulses as functions of time since "trigger" time. An analysis region is added to compute the total numbers of photons inside it.