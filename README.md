<p align="center"> <img src="./plots/s3p3_logo.png" width="30%"> </p>

# Transmission and reflection of light in PTFE for Dark Matter research

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17200120.svg)](https://doi.org/10.5281/zenodo.17200120)
[![GitHub license](https://img.shields.io/github/license/simop07/swiss_summer_student_program)](https://github.com/simop07/swiss_summer_student_program/blob/main/LICENSE)

## Overview
This repository contains material developed during my time as a Visiting Student at the [University of Zurich](https://www.uzh.ch/en.html) (UZH) as part of the [Swiss Summer Student Particle Physics Program](https://swiss.sspp.program.phys.ethz.ch/), held from 14 July to 14 September 2025.

The work from this project has been fully described in:

- [Presentation](https://indico.cern.ch/event/1579115/#30-transmission-and-reflection) $\rightarrow$ Overview of methodology, measurements, and key results.
- [Report](PasquiniSimoneReportSSSPPP.pdf) $\rightarrow$ Detailed documentation of experimental setup, analysis methods, and results.

## Summary
Polytetrafluoroethylene (PTFE), commonly known by its brand name Teflon, is extensively employed in rare event searches, including Dark Matter (DM) research. Indeed, in DM experiments, PTFE reflectors are widely used to surround Liquid Xenon (LXe) detectors, as PTFE excellent diffuse reflectivity allows to efficiently collect the Vacuum Ultra-Violet scintillation light ($\sim178$ nm) emitted by LXe. The critical influence of reflectivity on detector performance, along with the growing interest in assembling thinner PTFE detector walls that maintain high collection efficiency while minimizing light leakages, motivates the need for accurate knowledge of PTFE transmittance $(T)$ and reflectance $(R)$, and of their dependence on PTFE thickness. In this work, a `C++`/`ROOT` pulse analysis software is developed to determine $T$ and $R$ at the PTFE--air interface. Measurements are performed by using a blue LED light source and PTFE discs of varying thicknesses (from $0.20\pm0.05$ mm to $5.15\pm0.05$ mm) at incidence and reflection angles of $30^\circ$, $45^\circ$ and $60^\circ$. Attenuation coefficients are consistent across angles ($\lambda_T^{30^\circ}=(1.32\pm0.12)$ mm, $\lambda_T^{45^\circ}=(1.33\pm0.12)$ mm and $\lambda_T^{60^\circ}=(1.33\pm0.13)$ mm), while $R$ shows a strong angular dependence. These results will contribute to the future development of ALPINE detector, whose LXe TPC features a novel PTFE--metal--PTFE interface with light-tight properties yet to be comprehensively characterized.

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