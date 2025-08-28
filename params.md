# Params info

## Primary digitised data for $45^{\circ}$
### $1$ PTFE layer
| Region     | Inc_T [PE/ns] | Inc_R [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|---------------|---------------|----------------|--------------|--------|--------|
| PreTrig    | 0.000         | 0.000         | 0.000          | 0.000        | -nan   | -nan   |
| Trig       | 124.047       | 28.114        | 26.600         | 21.339       | 0.214  | 0.172  |
| PostTrig1  | 28.011        | 7.465         | 7.117          | 6.520        | 0.254  | 0.233  |
| PostTrig2  | 24.283        | 5.056         | 6.338          | 4.176        | 0.261  | 0.172  |

### $2$ PTFE layer
| Region     | Inc_T [PE/ns] | Inc_R [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|---------------|---------------|----------------|--------------|--------|--------|
| PreTrig    | 0.000         | 0.000         | 0.000          | 0.000        | -nan   | -nan   |
| Trig       | 124.047       | 28.114        | 15.240         | 25.203       | 0.123  | 0.203  |
| PostTrig1  | 28.011        | 7.465         | 4.719          | 6.952        | 0.168  | 0.248  |
| PostTrig2  | 24.283        | 5.056         | 4.084          | 6.530        | 0.168  | 0.269  |

### $3$ PTFE layer
| Region     | Inc_T [PE/ns] | Inc_R [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|---------------|---------------|----------------|--------------|--------|--------|
| PreTrig    | 0.000         | 0.000         | 0.000          | 0.000        | -nan   | -nan   |
| Trig       | 124.047       | 28.114        | 7.968          | 26.945       | 0.064  | 0.217  |
| PostTrig1  | 28.011        | 7.465         | 2.470          | 3.353        | 0.088  | 0.120  |
| PostTrig2  | 24.283        | 5.056         | 1.745          | 3.749        | 0.072  | 0.154  |

### Comparison table
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 1.55           | 0.214  | 0.172  |
| 3.10           | 0.123  | 0.203  |
| 5.14           | 0.064  | 0.217  |

## Secondary digitised data for $45^{\circ}$ (more points)
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.519  | 0.112  |
| 0.800          | 0.228  | 0.157  |
| 1.550          | 0.209  | 0.209  |
| 2.050          | 0.145  | 0.229  |
| 3.100          | 0.123  | 0.161  |
| 3.600          | 0.102  | 0.234  |
| 4.100          | 0.067  | 0.222  |
| 5.150          | 0.074  | 0.209  |

## Tertiary digitised data for $45^{\circ}$ (more points, more waveforms and with corrected geometry, with final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.506  | 0.143  |
| 0.800          | 0.236  | 0.292  |
| 1.550          | 0.210  | 0.426  |
| 2.050          | 0.143  | 0.427  |
| 3.100          | 0.115  | 0.520  |
| 3.600          | 0.095  | 0.511  |
| 4.100          | 0.069  | 0.567  |
| 5.150          | 0.071  | 0.572  |
**Fit results**  
- $\lambda_t ~ 0.7/0.8

## Quaternary digitised data for $30^{\circ}$  
*(fit transmittance up to 4.2 mm, 30k waveforms, corrected geometry with $\pm25^{\circ}$ acceptance, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.487  | 0.172  |
| 0.800          | 0.237  | 0.384  |
| 1.550          | 0.208  | 0.474  |
| 2.050          | 0.130  | 0.432  |
| 3.100          | 0.105  | 0.456  |
| 3.600          | 0.083  | 0.468  |
| 4.100          | 0.073  | 0.449  |
| 5.150          | 0.060  | 0.428  |
**Fit results**  
Chi2        = 0.004  
NDf         = 4  
Edm         = 0.000  
NCalls      = 79  
A           = 0.496 ± 0.052  
$\lambda_t$ = 0.849 ± 0.230  
B           = 0.084 ± 0.024  

## Quaternary digitised data for $45^{\circ}$  
*(fit transmittance up to 4.2 mm, 30k waveforms, corrected geometry with $\pm25^{\circ}$ acceptance, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.592  | 0.239  |
| 0.800          | 0.287  | 0.315  |
| 1.550          | 0.241  | 0.505  |
| 2.050          | 0.166  | 0.422  |
| 3.100          | 0.121  | 0.443  |
| 3.600          | 0.100  | 0.456  |
| 4.100          | 0.079  | 0.485  |
| 5.150          | 0.068  | 0.464  |
**Fit results**  
Chi2        = 0.005  
NDf         = 4  
Edm         = 0.000  
NCalls      = 72  
A           = 0.609 ± 0.060  
$\lambda_t$ = 0.842 ± 0.214  
B           = 0.098 ± 0.028  

## Quaternary digitised data for $60^{\circ}$  
*(fit transmittance up to 4.2 mm, 30k waveforms, corrected geometry with $\pm25^{\circ}$ acceptance, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 1.462  | 0.531  |
| 0.800          | 0.660  | 1.009  |
| 1.550          | 0.599  | 1.470  |
| 2.050          | 0.371  | 0.830  |
| 3.100          | 0.296  | 1.083  |
| 3.600          | 0.221  | 0.527  |
| 4.100          | 0.182  | 1.467  |
| 5.150          | 0.161  | 1.150  |
**Fit results**  
Chi2        = 0.047  
NDf         = 4  
Edm         = 0.000  
NCalls      = 77  
A           = 1.531 ± 0.187  
$\lambda_t$ = 0.795 ± 0.249  
B           = 0.236 ± 0.080  

## 5th data taking
*(fit transmittance up to 4.2 mm, 60k waveforms, corrected geometry with $\pm25^{\circ}$ acceptance, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

## Results in trigger region for 30°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.471  | 0.212  |
| 0.800          | 0.221  | 0.364  |
| 1.550          | 0.196  | 0.394  |
| 2.050          | 0.132  | 0.417  |
| 3.100          | 0.103  | 0.485  |
| 3.600          | 0.080  | 0.393  |
| 4.100          | 0.067  | 0.487  |
| 5.150          | 0.056  | 0.415  |
**Fit results for 30°**  
- $A = 0.481$  
- $\lambda_t = 0.811$  
- $B = 0.083$  


## Results in trigger region for 45°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.603  | 0.189  |
| 0.800          | 0.279  | 0.396  |
| 1.550          | 0.235  | 0.342  |
| 2.050          | 0.169  | 0.349  |
| 3.100          | 0.121  | 0.369  |
| 3.600          | 0.099  | 0.414  |
| 4.100          | 0.080  | 0.416  |
| 5.150          | 0.069  | 0.403  |
**Fit results for 45°**  
- $A = 0.629$  
- $\lambda_t = 0.779$  
- $B = 0.103$  


## Results in trigger region for 60°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 1.448  | 0.430  |
| 0.800          | 0.619  | 0.848  |
| 1.550          | 0.578  | 1.311  |
| 2.050          | 0.367  | 0.783  |
| 3.100          | 0.288  | 0.995  |
| 3.600          | 0.213  | 0.569  |
| 4.100          | 0.174  | 1.536  |
| 5.150          | 0.146  | 1.038  |
**Fit results for 60°**  
- $A = 1.544$  
- $\lambda_t = 0.732$  
- $B = 0.239$  

## 6th data taking
*(fit transmittance up to 4.2 mm, 60k waveforms, corrected geometry with $\pm25^{\circ}$ acceptance and correction with C(x) function, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

## Results in trigger region for 30°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.471  | 0.255  |
| 0.800          | 0.221  | 0.437  |
| 1.550          | 0.196  | 0.474  |
| 2.050          | 0.132  | 0.500  |
| 3.100          | 0.103  | 0.583  |
| 3.600          | 0.080  | 0.472  |
| 4.100          | 0.067  | 0.585  |
| 5.150          | 0.056  | 0.499  |

## Results in trigger region for 45°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.603  | 0.218  |
| 0.800          | 0.279  | 0.456  |
| 1.550          | 0.235  | 0.394  |
| 2.050          | 0.169  | 0.402  |
| 3.100          | 0.121  | 0.425  |
| 3.600          | 0.099  | 0.477  |
| 4.100          | 0.080  | 0.479  |
| 5.150          | 0.069  | 0.464  |

## Results in trigger region for 60°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 1.448  | 0.473  |
| 0.800          | 0.619  | 0.932  |
| 1.550          | 0.578  | 1.442  |
| 2.050          | 0.367  | 0.861  |
| 3.100          | 0.288  | 1.094  |
| 3.600          | 0.213  | 0.626  |
| 4.100          | 0.174  | 1.690  |
| 5.150          | 0.146  | 1.141  |
**Fit results for 30°**
-A           = 0.481 ± 0.055
-lambda_t    = 0.811 ± 0.239 mm
-B           = 0.083 ± 0.024

**Fit results for 45°**
-A           = 0.629 ± 0.067
-lambda_t    = 0.779 ± 0.211 mm
-B           = 0.103 ± 0.028

**Fit results for 60°**
-A           = 1.544 ± 0.212
-lambda_t    = 0.732 ± 0.251 mm
-B           = 0.239 ± 0.080

-Correction factor for 30° = 0.462
-Correction factor for 45° = 0.472
-Correction factor for 60° = 0.505

-Correction factor with geo correction for 30° = 0.385
-Correction factor with geo correction for 45° = 0.410
-Correction factor with geo correction for 60° = 0.459

-Fit result 30°: p0 = 0.641 ± 0.022
-Fit result 45°: p0 = 0.621 ± 0.037
-Fit result 60°: p0 = 1.512 ± 0.143

## 7th data taking
*(fit transmittance up to 5.3 mm, 60k waveforms, corrected geometry with $\pm25^{\circ}$ acceptance and correction with C(x) function, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

### Results in trigger region for 30°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.471  | 0.255  |
| 0.800          | 0.221  | 0.437  |
| 1.550          | 0.196  | 0.474  |
| 2.050          | 0.132  | 0.500  |
| 3.100          | 0.103  | 0.583  |
| 3.600          | 0.080  | 0.472  |
| 4.100          | 0.067  | 0.585  |
| 5.150          | 0.056  | 0.499  |

### Results in trigger region for 45°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 0.603  | 0.218  |
| 0.800          | 0.279  | 0.456  |
| 1.550          | 0.235  | 0.394  |
| 2.050          | 0.169  | 0.402  |
| 3.100          | 0.121  | 0.425  |
| 3.600          | 0.099  | 0.477  |
| 4.100          | 0.080  | 0.479  |
| 5.150          | 0.069  | 0.464  |

### Results in trigger region for 60°  
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 0.200          | 1.448  | 0.473  |
| 0.800          | 0.619  | 0.932  |
| 1.550          | 0.578  | 1.442  |
| 2.050          | 0.367  | 0.861  |
| 3.100          | 0.288  | 1.094  |
| 3.600          | 0.213  | 0.626  |
| 4.100          | 0.174  | 1.690  |
| 5.150          | 0.146  | 1.141  |
**Fit results for 30°**
-A           = 0.480 ± 0.049
-lambda_t    = 0.891 ± 0.226 mm
-B           = 0.072 ± 0.020

**Fit results for 45°**
-A           = 0.627 ± 0.061
-lambda_t    = 0.849 ± 0.201 mm
-B           = 0.091 ± 0.023

**Fit results for 60°**
-A           = 1.533 ± 0.186
-lambda_t    = 0.813 ± 0.241 mm
-B           = 0.206 ± 0.068

-Correction factor for 30° = 0.462
-Correction factor for 45° = 0.472
-Correction factor for 60° = 0.505

-Correction factor with geo correction for 30° = 0.385
-Correction factor with geo correction for 45° = 0.410
-Correction factor with geo correction for 60° = 0.459

-Fit result 30°: p0 = 0.641 ± 0.022
-Fit result 45°: p0 = 0.621 ± 0.037
-Fit result 60°: p0 = 1.512 ± 0.143

## 8th data taking
*(fit transmittance up to 5.3 mm, 60k waveforms, 2 sigmas in the trigger region, corrected geometry with $\pm25^{\circ}$ acceptance and correction with C(x) function, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

### 30°  
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200          | 0.469  | 0.004| 0.227  | 0.006|
| 0.800          | 0.223  | 0.003| 0.429  | 0.008|
| 1.550          | 0.198  | 0.002| 0.476  | 0.008|
| 2.050          | 0.131  | 0.002| 0.504  | 0.008|
| 3.100          | 0.105  | 0.002| 0.570  | 0.009|
| 3.600          | 0.080  | 0.001| 0.547  | 0.009|
| 4.100          | 0.065  | 0.001| 0.584  | 0.009|
| 5.150          | 0.057  | 0.001| 0.592  | 0.009|

### 45°  
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200          | 0.613  | 0.006| 0.219  | 0.012|
| 0.800          | 0.283  | 0.004| 0.404  | 0.014|
| 1.550          | 0.250  | 0.003| 0.530  | 0.015|
| 2.050          | 0.165  | 0.002| 0.541  | 0.015|
| 3.100          | 0.124  | 0.002| 0.448  | 0.014|
| 3.600          | 0.100  | 0.002| 0.496  | 0.014|
| 4.100          | 0.083  | 0.002| 0.628  | 0.015|
| 5.150          | 0.071  | 0.001| 0.492  | 0.014|

### 60°  
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200          | 1.406  | 0.022| 0.503  | 0.062|
| 0.800          | 0.613  | 0.011| 0.396  | 0.060|
| 1.550          | 0.559  | 0.010| 0.859  | 0.065|
| 2.050          | 0.346  | 0.008| 0.806  | 0.065|
| 3.100          | 0.288  | 0.006| 0.933  | 0.067|
| 3.600          | 0.217  | 0.005| 0.565  | 0.063|
| 4.100          | 0.173  | 0.005| 1.499  | 0.073|
| 5.150          | 0.151  | 0.005| 0.899  | 0.067|
**Fit results**
### 30°  
- A           = 0.433 ± 0.005  
- $\lambda_t$ = 1.213 ± 0.022 mm  
- B           = 0.055 ± 0.001  

### 45°  
- A           = 0.552 ± 0.007  
- $\lambda_t$ = 1.186 ± 0.025 mm  
- B           = 0.069 ± 0.002  

### 60°  
- A           = 1.220 ± 0.025  
- $\lambda_t$ = 1.162 ± 0.040 mm  
- B           = 0.153 ± 0.005  

### Correction Factors
| Angle | Correction | Correction with Geo |
|-------|------------|-------------------|
| 30°   | 0.462      | 0.385             |
| 45°   | 0.472      | 0.410             |
| 60°   | 0.505      | 0.459             |

### Additional Fit Results
| Angle | p0 |
|-------|------|
| 30°   | 0.657 ± 0.004 |
| 45°   | 0.671 ± 0.006 |
| 60°   | 1.239 ± 0.026 |

## 9th data taking
*(fit transmittance up to 5.3 mm, 60k waveforms, 1.8 sigmas in the trigger region, corrected geometry with $\pm25^{\circ}$ acceptance and correction with C(x) function, final update on [waveformAnalysis.cpp](waveformAnalysis.cpp) parameters)*

### Results in trigger region for 30°
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200 | 0.468 | 0.004 | 0.239 | 0.013 |
| 0.800 | 0.221 | 0.003 | 0.433 | 0.020 |
| 1.550 | 0.196 | 0.002 | 0.479 | 0.022 |
| 2.050 | 0.131 | 0.002 | 0.507 | 0.023 |
| 3.100 | 0.104 | 0.002 | 0.576 | 0.025 |
| 3.600 | 0.079 | 0.001 | 0.553 | 0.024 |
| 4.100 | 0.065 | 0.001 | 0.589 | 0.025 |
| 5.150 | 0.057 | 0.001 | 0.608 | 0.026 |

### Results in trigger region for 45°
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200 | 0.613 | 0.006 | 0.223 | 0.017 |
| 0.800 | 0.282 | 0.004 | 0.421 | 0.022 |
| 1.550 | 0.251 | 0.003 | 0.514 | 0.025 |
| 2.050 | 0.164 | 0.002 | 0.532 | 0.026 |
| 3.100 | 0.124 | 0.002 | 0.454 | 0.023 |
| 3.600 | 0.100 | 0.002 | 0.504 | 0.025 |
| 4.100 | 0.083 | 0.002 | 0.617 | 0.028 |
| 5.150 | 0.070 | 0.001 | 0.496 | 0.024 |

### Results in trigger region for 60°
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200 | 1.406 | 0.022 | 0.484 | 0.070 |
| 0.800 | 0.616 | 0.011 | 0.333 | 0.066 |
| 1.550 | 0.560 | 0.010 | 0.746 | 0.075 |
| 2.050 | 0.346 | 0.008 | 0.804 | 0.077 |
| 3.100 | 0.288 | 0.006 | 0.871 | 0.079 |
| 3.600 | 0.219 | 0.005 | 0.576 | 0.072 |
| 4.100 | 0.171 | 0.005 | 1.447 | 0.092 |
| 5.150 | 0.151 | 0.005 | 0.851 | 0.079 |

### Results in trigger region for 61°
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200 | 1.438 | 0.023 | 0.273 | 0.066 |
| 0.800 | 0.630 | 0.012 | 0.686 | 0.077 |
| 1.550 | 0.527 | 0.010 | 0.740 | 0.078 |
| 2.050 | 0.359 | 0.008 | 0.586 | 0.073 |
| 3.100 | 0.290 | 0.008 | 0.843 | 0.080 |
| 3.600 | 0.237 | 0.007 | 0.599 | 0.074 |
| 4.100 | 0.179 | 0.005 | 1.182 | 0.089 |
| 5.150 | 0.150 | 0.004 | 0.988 | 0.084 |

### Results in trigger region for 62°
| Thickness [mm] | Prob_T | ErrT | Prob_R | ErrR |
|----------------|--------|------|--------|------|
| 0.200 | 1.486 | 0.022 | 0.530 | 0.070 |
| 0.800 | 0.668 | 0.011 | 0.757 | 0.066 |
| 1.550 | 0.546 | 0.010 | 0.902 | 0.075 |
| 2.050 | 0.380 | 0.008 | 1.031 | 0.077 |
| 3.100 | 0.303 | 0.006 | 1.059 | 0.079 |
| 3.600 | 0.223 | 0.005 | 0.914 | 0.072 |
| 4.100 | 0.183 | 0.005 | 0.971 | 0.092 |
| 5.150 | 0.146 | 0.005 | 0.948 | 0.079 |

---

### Fit results
**30°**
- A = 0.432 ± 0.005  
- λₜ = 1.203 ± 0.022 mm  
- B = 0.055 ± 0.001  

**45°**
- A = 0.551 ± 0.007  
- λₜ = 1.190 ± 0.026 mm  
- B = 0.069 ± 0.002  

**60°**
- A = 1.224 ± 0.026  
- λₜ = 1.172 ± 0.040 mm  
- B = 0.151 ± 0.006  

**61°**
- A = 1.265 ± 0.028  
- λₜ = 1.136 ± 0.038 mm  
- B = 0.153 ± 0.005  

**62°**
- A = 1.350 ± 0.026  
- λₜ = 1.210 ± 0.040 mm  
- B = 0.142 ± 0.006  

---

### Correction factors
- 30° = 0.462  
- 45° = 0.472  
- 60° = 0.505  
- 61° = 0.505  
- 62° = 0.505  

**With geo correction**
- 30° = 0.385  
- 45° = 0.410  
- 60° = 0.459  
- 61° = 0.459  
- 62° = 0.459  

---

### p0 Fit results
- 30°: p0 = 0.668 ± 0.008  
- 45°: p0 = 0.683 ± 0.009  
- 60°: p0 = 1.192 ± 0.030  
- 61°: p0 = 1.195 ± 0.031  
- 62°: p0 = 1.358 ± 0.020  

## Oscilloscope data for [4LayersA.txt](./data/4LayersA.txt)
```
  // Bin for hPhotoelectron = 60

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.7, 1.5);

  // Make the average of first 1100 samples of a waveform
  void baseline(int nInitialSamples = 1100);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 1100);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 55 mV*ns)
  double const areaConvFactor{55.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{1.1410014154 - 2 * 0.2461653408};
  double const triggerEnd{1.1410014154 + 2 * 0.2461653408};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.1};
  double const preTriggerEnd{0.7};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{1.7};
  double const postTriggerEnd1{2.3};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{3.};
  double const postTriggerEnd2{3.6};
```

## Oscilloscope data for $30^{\circ}$
| Region     | Inc [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|-------------|----------------|--------------|--------|--------|
| PreTrig    | 1.897       | 1.700          | 2.387        | 0.896  | 1.258  |
| Trig       | 1055.665    | 229.932        | 380.482      | 0.218  | 0.360  |
| PostTrig1  | 472.744     | 117.550        | 150.555      | 0.249  | 0.318  |
| PostTrig2  | 254.758     | 84.519         | 74.242       | 0.332  | 0.291  |

### [incT.txt](./data/incT.txt)
```
  // Bin for hPhotoelectron = 60

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);
  
  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA0.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{37.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5742169474 - 2 * 0.0945427226};
  double const triggerEnd{0.5742169474 + 2 * 0.0945427226};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.45};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.8};
  double const postTriggerEnd1{1.2};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.4};
  double const postTriggerEnd2{1.8};
```

### [transm.txt](./data/transm.txt)
```
  // Bin for hPhotoelectron = 60

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA2.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{37.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5880322313 - 2 * 0.0927908366};
  double const triggerEnd{0.5880322313 + 2 * 0.0927908366};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.5};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.85};
  double const postTriggerEnd1{1.25};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.45};
  double const postTriggerEnd2{1.85};
```

### [refl.txt](./data/refl.txt)
```
  // Bin for hPhotoelectron = 25

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA3.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{60.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5745826154 - 2 * 0.1743512132};
  double const triggerEnd{0.5745826154 + 2 * 0.1743512132};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.5};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.95};
  double const postTriggerEnd1{1.35};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.45};
  double const postTriggerEnd2{1.85};
```

## Oscilloscope data for $60^{\circ}$ ~ WORK IN PROGRESS
| Region     | Inc [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|-------------|----------------|--------------|--------|--------|
| PreTrig    | 1.147       | 0.000          | 3.370        | 0.000  | 2.938  |
| Trig       | 243.475     | 215.317        | 1043.009     | 0.884  | 4.284  |
| PostTrig1  | 174.684     | 101.659        | 439.060      | 0.582  | 2.513  |
| PostTrig2  | 111.794     | 80.781         | 189.913      | 0.723  | 1.699  |

### [inc60T.txt](./data/inc60T.txt)
```
  // Info hPhotoelectron
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               40, 0, 4.5);

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA2.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{27.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5670213385 - 2 * 0.0938516294};
  double const triggerEnd{0.5670213385 + 2 * 0.0938516294};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.45};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.8};
  double const postTriggerEnd1{1.2};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.4};
  double const postTriggerEnd2{1.8};
```

### [transm60.txt](./data/transm60.txt)
```
  // Info hPhotoelectron
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               40, 0, 4.5);

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA2.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{27.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5612464817 - 2 * 0.0749055716};
  double const triggerEnd{0.5612464817 + 2 * 0.0749055716};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.5};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.85};
  double const postTriggerEnd1{1.25};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.45};
  double const postTriggerEnd2{1.85};
```

### [refl60.txt](./data/refl60.txt)
```
  // Info for hPhotoelectron
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               40, 0, 4.5);

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);
  
  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA60_2.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{70.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5758587068 - 2 * 0.0889747690};
  double const triggerEnd{0.5758587068 + 2 * 0.0889747690};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.45};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.8};
  double const postTriggerEnd1{1.2};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.4};
  double const postTriggerEnd2{1.8};
```