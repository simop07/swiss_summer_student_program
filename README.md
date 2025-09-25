<p align="center"> <img src="./plots/s3p3_logo.png" width="30%"> </p>

# UZH Swiss Summer Student Particle Physics Program

[![DOI](https://zenodo.org/badge/1024295677.svg)](https://doi.org/10.5281/zenodo.16323819)
[![GitHub license](https://img.shields.io/github/license/simop07/swiss_summer_student_program)](https://github.com/simop07/swiss_summer_student_program/blob/main/LICENSE)

## Overview
This repository contains material developed during my time as a Visiting Student at the [University of Zurich](https://www.uzh.ch/en.html) (UZH) as part of the [Swiss Summer Student Program](https://swiss.sspp.program.phys.ethz.ch/) in Particle Physics, held from 14 July to 14 September 2025.

The work from this project has been fully described in:

- [Presentation](https://indico.cern.ch/event/1579115/#30-transmission-and-reflection) $\rightarrow$ Overview of methodology, measurements, and key results.
- [Report](PasquiniSimoneReportSSSPPP.pdf) $\rightarrow$ Detailed documentation of experimental setup, analysis methods, and results.

## Summary
Polytetrafluoroethylene (PTFE), often known by its brand name Teflon, is extensively used in rare event searches, including Dark Matter (DM) experiments. Indeed, to efficiently collect the Vacuum Ultra-Violet scintillation light ($\sim178$ nm) emitted by Liquid Xenon (LXe) detectors employed in DM research, PTFE reflectors are commonly utilized to surround LXe volumes due to their excellent diffuse reflectivity. The critical influence of reflectivity on detector performance, along with the growing interest in assembling thinner PTFE detector walls that maintain high collection efficiency while minimizing light leakages, motivates the need for accurate knowledge of PTFE transmittance $(T)$ and reflectance $(R)$ dependence on PTFE thickness. In this work, a `C++`/`ROOT` pulse analysis software is developed to determine $T$ and $R$ at the PTFE--air interface, by using a blue LED light source and PTFE discs of varying thicknesses (from $0.2\pm0.05$ to $5.15\pm0.05$ mm). Measurements are performed at incidence and reflection angles of $30^\circ$, $45^\circ$ and $60^\circ$. Attenuation coefficients are consistent across angles ($\lambda_T^{30^\circ}=(1.32\pm0.12)$ mm, $\lambda_T^{45^\circ}=(1.33\pm0.12)$ mm and $\lambda_T^{60^\circ}=(1.33\pm0.13)$ mm), while $R$ shows a strong angular dependence. These results will contribute to the future development of ALPINE detector, whose LXe TPC features a novel PTFE--metal--PTFE interface with light-tight properties yet to be comprehensively characterized.

## Pulse analysis
The pulse analysis software extracts key parameters from photomultiplier (PMT) signals to determine PTFE transmittance and reflectance. The workflow includes:

- Processing digitized PMT waveforms.
- Identifying and integrating pulses.
- Applying geometric acceptance corrections.
- Computing transmittance and reflectance for different PTFE thicknesses and angles.

## Repository Contents
- [`merge.cpp`](src/merge.cpp) $\rightarrow$ Merge multiple waveform files into a single dataset.
- [`compare.cpp`](src/compare.cpp) $\rightarrow$ Compare 1D histograms for the analysis.
- [`lightAnalysis.cpp`](src/lightAnalysis.cpp) $\rightarrow$ Compute PTFE transmittance and reflectance.
- [`convertCSVtoTXT.cpp`](src/convertCSVtoTXT.cpp) $\rightarrow$ Convert files from .CSV to .txt.
- [`waveformAnalysisPos.cpp`](src/waveformAnalysisPos.cpp) / [`waveformAnalysisNeg.cpp`](src/waveformAnalysisNeg.cpp) $\rightarrow$ Custom waveform analysis devoted to analyze pulse estraction and analysis.
- [`waveform.cpp`](src/waveform.cpp) / [`waveformOsc.cpp`](src/waveformOsc.cpp) / [`waveformOscDownsample.cpp`](src/waveformOscDownsample.cpp) $\rightarrow$ Process digitized or oscilloscope data.