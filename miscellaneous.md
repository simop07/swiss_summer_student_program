# Miscellaneous notes

## 15/07/2025
We want to send a pulse signal to an **LED** and observe the response from the photomultiplier tube (PMT). However, the pulse signal isn't just a simple amplitude - it consists of two parts: a **baseline** (a logical ‘0' which is not at $0$ $V$, set at $3$ $V$) and a peak (fixed at $4.5$ $V$). What we're trying to do is minimize the difference between the peak and the baseline, to minimize noise in the system.

### System Overview
- The LED is driven by a pulse signal, meaning it's only turned on during short bursts (the peaks).
- Light from the LED passes through a PTFE sample (used as an attenuator/diffuser). Probably we'll refer to this as "filters" in the following.
- The PMT detects the light and produces an analog electrical signal – of course it is connected to a HV in order to amplify the light yield.
- We want to measure this PMT output to understand how much light is reaching it.
- The pulse signal that drives the LED is used both to control the LED and as a trigger for acquisition - so we're sure the PMT's signal is acquired only when light is being emitted. Note that the trigger we're using is the LED signal transformed into TTL logic.

### Pulse Signal Shape
- Baseline: the resting voltage level when the LED is off.
- Peak: the voltage level during the pulse - it turns the LED on.
- The difference between peak and baseline (called pulse height or amplitude) determines how much current goes to the LED, and therefore how much light is emitted.

### What Are We Trying to Optimize?
Our goal is to minimize the noise in the PMT's signal. More specifically:
- Keep the baseline flat and noise-free.
- Have a peak just strong enough to make the LED emit a detectable, measurable amount of light.
- Avoid making the peak too strong, which would cause:
- The LED to emit too much light (we prefer it to work above threshold), but at the same time we don't want to use all the light that it emits to prevent:
    - The PMT to go into saturation
    - The system to become non-linear

So, it's not about minimizing the pulse height to zero - it's about **finding the smallest pulse height that gives a clear, stable signal well above the noise, without overdriving the LED or PMT**.

### What Noise Are We Talking About?
The noise here refers mostly to:
- Electronic noise in the PMT and the acquisition electronics
- Baseline fluctuations when the LED is off
- Background light or reflections that may still trigger the PMT

If the LED pulse is too strong, we may:
- Introduce electrical ringing
- Increase afterpulses in the PMT
- Degrade the signal-to-noise ratio (SNR) due to PMT saturation
Thus, a larger pulse implies stronger LED light, which in turn implies possibly higher background noise and non-linear response.

We sort of want that the PMT collects all the possible light, so we choose a proper value of HV in order to collect all the light produced by the LED. So this seems a sort of "driven saturation", but it is actually different from the saturation mentioned earlier. Indeed, even if HV is fixed, if the LED emits too much light, the PMT can become saturated because:
- Too many photoelectrons are generated too quickly
- The multiplication chain (dynodes) gets overwhelmed
- The output current becomes non-linear (i.e. no longer proportional to the number of incident photons)
This leads to:
- Signal flattening or clipping (looks like a "flat top" on the oscilloscope)
- Loss of linearity (you can't trust the signal amplitude anymore)
- Pulse pile-up or distortion

# 16/07/2025
## Dark Matter: Historical Evidence and Detection
### Zwicky and the Virial Theorem
- Using thermodynamics and the virial theorem, Zwicky analyzed galaxy dynamics: by observing how things move, we can estimate how massive a galaxy must be to remain stable.
- Observation: In the outer regions of galaxies, stars move faster than expected based on visible mass.

### Vera Rubin (1960s)
- Rubin emphasized that this anomaly was not limited to a single galaxy. Rotation curves across many galaxies showed stars moving too fast in the outskirts, implying the presence of unseen mass.
- She framed the problem on a cosmological scale, arguing that much more mass must exist than what we observe.

## Initial Interpretations of Dark Matter (DM):
- Maybe we don't fully understand gravity (e.g. MOND, Modified Newtonian Dynamics – not widely supported due to fine-tuning issues).
- Or: Dark Matter exists.

## The Bullet Cluster (Coma Cluster)
A powerful observational proof: collision of two galaxy clusters.
- Gravitational lensing reveals the mass distribution (blue areas), showing where gravity "acts".
- Hot X-ray-emitting gas (visible matter) lags behind. This spatial separation shows most of the mass isn't in the gas, but rather in something invisible and collisionless: Dark Matter.

## Properties Dark Matter must have:
- Neutral (no electromagnetic interaction)
- Massive (contributes to gravitational potential)
- Cold (non-relativistic; otherwise structure formation would be suppressed – e.g., sterile neutrinos would wash out galactic structures)
- Stable over the lifetime of the Universe (still abundant today)

We know the total amount of DM from cosmological measurements (CMB, structure formation), and estimating its mass is key to quantifying how much of it is present.

## Radiation: concepts and detection
- Radiation = Energy Transfer between particles.
- We distinguish between:
	- Ionising radiation (can remove electrons from atoms): alpha, beta, neutrons, neutrinos
	- Non-ionising radiation (e.g. radio, microwave)
- Focus on ionising types (used in DM detection):
	- Alpha particles, beta particles, neutrons, neutrinos

## Dark Matter detection requirements
We aim to detect DM particles (non-relativistic, neutral) moving at ~400 km/s.
- Photons and beta particles are ultra-relativistic
- Alphas and neutrons are non-relativistic (neutrons often used as DM simulants because they're neutral)
- Typical neutron speeds: few hundred m/s
When a DM particle interacts with ordinary matter, it may produce:
- Heat
- Light
- Charge
These are not independent - they are connected through energy conversion processes. E.g., ionisation (charge) can lead to recombination (light), and motion can lead to vibrations which further leads to heat.

## Triangle of radiation detection
1. Light Detection
	- Scintillation: energy deposit implying excited states implying ion-electron pairs implying recombination implying photons
	- Many scintillators are not transparent to their own scintillation light; hence, doping (impurities) is added
	- Light sensors include PMTs, SiPMs, etc.
2. Charge Detection
	- Proportional counters: Gas-filled tubes with anode wire at high voltage
	- Ionization implying electrons drift to anode; positive ions to cathode
	- Semiconductor detectors (e.g. CCD, germanium):
	- CCD reads position domain, silicon diodes read in time domain
3. Heat Detection
	- Temperature changes are measured as energy deposition
	- Requires operation at millikelvin temperatures (10 mK) implying thermal baths
	- Individual particle interactions cause measurable temperature spikes

## Types of interactions
- Gamma rays can scatter off electrons or nuclei, producing ionisation and heat
- Important to distinguish electron recoils from nuclear recoils
Measuring multiple observables simultaneously improves signal discrimination (already done in nuclear physics using scintillators).

## The XENON Experiment (Dark Matter Direct Detection)
- Goal: Large, sensitive detector with high target mass implying Liquid Xenon (LXe) is ideal
	- High density, scalable, purifiable
	- Self-shielding
	- No isotopes, so minimal background
	- No dopants needed: excited LXe forms excimers that emit VUV scintillation light
	- LXe is transparent to its own scintillation (property of noble liquids), and this happens because the light is not emitted directly by Xenon nuclei, but by Xenon excimers! Emitted light is then not absorbed by ground state of Xenon nuclei.
- Detector Structure:
	- 7 tons of LXe
	- 500 light sensors (PMTs): half on top, half at the bottom
	- Two-phase detector: liquid + gas
	- Electric field applied via metal wires
- Signal Types:
	- S1: Prompt scintillation in liquid (photons)
	- S2: Delayed signal from ionised electrons that drift into the gas and generate secondary scintillation light
- Each electron produces 30ish photons in S2 implying high sensitivity
- S2/S1 ratio is larger for electron recoils than for nuclear recoils implying helps identify event type
- Time difference between S1 and S2 implying depth of interaction
- Hotspot location on sensor array implying planar (x-y) position implying position resolution down to a few mm
- Fiducialization: Using reconstructed position to select only events from the central volume (lowest background)
- Neutrinos: can also scatter off nuclei implying potential background for DM detection
- Outer Layers:
	- LXe Skin: outer volume used as a veto, looking for coincident signals that may indicate background particles entering the detector

# 17/07/2025
- Aim: measure reflectance and absorbance of light in PTFEs by observing the scintillation signal size detected by PMTs when shining light from a blue LED.
- Experimental setup:
    - Pulse generator used to light a blue LED with:
    - Frequency: 5000 Hz
    - High Level: 4.0 V
    - Low Level: 2.3 V
    - Delay: 0 ns
    - Width: 70.0 ns
    - Trail 3.0 ns
    - High voltage generator set at 1100 V
    - PMTs
    - LED and PTFE holders. To improve the geometry of the LED holder, a better option was to simply leave one side of the box open so that Teflon could be simply inserted in instead of "being pushed" into the cave box.
    - Wires to connect everything – note that the trigger happens on the TTL logic cable correspondent with the pulse signal.

Initially the idea is to detect single photon events – so that the signal we will analyse is narrow (a larger signal will simply be more difficult to interpret (is it simply a delayed pulse of the same PMT event or something due to the reflectance?)). To do so, we need to reduce the width of the pulse as much as possible, and we need to find the limiting condition maximising the Low Level and minimising the High Level. Why?
Because the total number of photons emitted per pulse by the LED is proportional to the area under the pulse:
$$(V_{high}-V_{low} )*Width$$
Now the question is: why the Low Level should be as high as possible? Operating close to threshold means that the LED only turns on briefly, and often emits just 0 or 1 photons per pulse. If the Low Level is too low, the LED is strongly reverse biased before the pulse  when the pulse hits, the LED has a larger forward bias and emits much more light. In this way the spikes seen in the PMT readout are more similar before and after the thick forest of the readout – even if there's still some discrepancy between before and after the LED pulse. So the best option is putting Low Level at 2.3 V because at 2.3 V the LED is barely on.

Instead, as to the high level, it should be not too low because otherwise the signal to noise ratio will be too low, but it shouldn't be too high because in that case we produce LOTS of photons instead of just one. So we should minimise the High Level.

# 29/07/2025
Some **useful insights on the analysis that we're going to do next**. The digitizer has two important parameters that need to be tuned: the voltage range that needs to be digitized, like from 0 to -14 mV, and the number of bits used to digitise it (16, 32 bits…), which gives a smoother or rougher curve depending on the needs. Plus, when we will use the double PMT to measure both reflectance and absorption at the same time, we will need two sets of PTFEs: one in front of the LED (light source), whose aim is just to obtain 1 PE or 2 PE peaks, and the one in front of the PMT, whose aim is study the probabilities of reflectance and transmission dependencies on its thickness. Question is, why should we use the first set of PTFE? Adding more layers increases the chances that light will:
   - be reflected multiple times within the layers, coming back to the LED and being absorbed
   - be diffused more uniformly, and with more stability
   - and the light directed towards the PMT will be more controlled  obtain clearer 1 or 2 PE peaks.
This effectively filters out stray light, which is a common source of background noise in PMT measurements.

Before using the digitiser, I need to do the following:
- In the oscilloscope data, remove the waves which peak at -14 mV, because probably their area information is wrong (they would have a larger peak, but the oscilloscope (with our setup) only arrives at -14 mV). Modify the measurement range of the oscilloscope in the future if needed.
- What is happening in the width histogram? Try to understand this by plotting the code which Nicolas has shared with you (most interesting parameters are width vs area and noise vs area)
- Create the code involving reflectance and transmittance with themselves and as a function of the thickness of PTFE (see notes).
- Fixing the thickness, plot what happens at different angles (so transmittance and reflectance as a function of the incident angle).

# 05/08/2025
## [waveformOsc.cpp](waveformOsc.cpp) miscellaneous information
At the moment I've made just one selection of data by removing waveform reaching "maximum" peak value ($-13.9333$ $mV$).

## [waveform.cpp](waveform.cpp) miscellaneous information
At the moment I've made - it is commented - just one selection of data by removing pulses whose area in PE is under a certain value.

## Development of the pulse analysis for both [waveformOsc.cpp](waveformOsc.cpp) and [waveform.cpp](waveform.cpp)
After analysing all possible parameters of the analysis (check the above files to see them more specifically), I need to find a way to separate the noise pulses from the physical ones. To do so, I sum all the pulses and fit a gaussian in the area which I would define "the trigger region". Then, I identify other 3 regions, a pre-trigger region and two post-trigger regions. The two post trigger regions are built in such a way that the first one coincides exactly with the decay of the "sum of pulses" graph, while the other corresponds with the last part of the pulses time.

Subsequently, I found several parameters for each region such as:
- Pulses region / total (note that a pulse belongs to a region if its peaks belong to that region)
- Total area (simply the sum of the areas of the pulses belonging to a region)
- Total number of PE (area scaled with the conversion factor into PE)
- Number of PE per pulse (total number of PE / number of pulses within a region)
- Rate (total number of PE / over lifetime)
Then, the idea is to further clean the graphs adding selections. Selections are used to better discriminate between physical signals and noise signals. The idea is that not only we search for pulses using the threshold cointained in the waveformAnalysis class, but also imposing some constraints on the area. Again, the idea is that a signal is considered physical if its area is above a certain value. The ratio underlying this first selection can be summarised in the following expression (note that the right hand side is usually expressed in ACD counts, while the left hand side is an area!):
$$\mu_{singlePE}-N\cdot\sigma_1^1\gg\mu_{noise}+M\cdot\sigma_{noise}$$
meaning that the right tail of the noise should be much lower than the left tail of the single PE event.

# 06/08/2025
Some random comments I received at the meeting today:
- To investigate undershots for pulses, it is better to look at wide pulses: in fact, wider pulses are simply more prone to have more undershots.
- A nice thing that we could do is making the PTFE used for reflectance measurements a little be **wet** (be careful, you're using high voltage!!) -> in this way, you can mimick the interface Xenon-PTFE in a better way more than using air-PTFE. When performing meaurements of reflectance, you should also be aware of the **critic angle** at the interface.
- I am doing what I am doing for ALPINE detector -> for its segmentation. Basically we want to know if and how many photons will pass from one segmented layer to another.

# 11/08/2025
When in [lightAnalysis.cpp](lightAnalysis.cpp) i mention "Rate correct" I mean that I am considering the rate of a given region with a reduced number of PE: indeed, I am subtracting the PE of the pre-trigger region in the following way:
$$N_{PE region}^{correct}=N_{PE region}-R_{bkg}*\Delta t_{region}$$
$$RATE^{correct}=N_{PE region}^{correct}/\Delta t_{region}$$
where $R_{bkg}=N_{bkg}/\Delta t_{bkg}$.

# 15/08/2025
I further developed what I did [here](#15082025) in order to improve measurements on reflectance. Indeed, reflected photons reaching PMT_R include also incident photons that directly travel from the LED to the PMT_R. So what I did was take the $N_{PE region}^{correct}$ and subtract then the rate of incident photons $RATE_{IncR​}$, in the following way:
$$rateCorrRefl_{R}=RATE^{correct}-RATE_{IncR​}$$

# 21/08/2025
Another issue is **geometric acceptance**. PMTs see only a fraction set by the solid angle they subtend at the emission point. To compute a possible correction function, I've used what is in this paper [here](https://repositorium.uni-muenster.de/document/miami/42495273-2cc7-4632-99af-accf074921e7/diss_levy.pdf). There's a function (Eq. $6.17$) that models the reflections off a PTFE layer of 5 mm thickness. In particular, I'm analysing the function by choosing a specific value of 45° angle of incidence between the light beam and the normal to the PTFE surface. The code then computes the geometrical correction factor by integrating the angular response function F(x) over the PMT acceptance window (e.g. 40°-50°) and normalizing to the full angular range (–90° to 90°). In practice, this ratio represents the fraction of reflected photons effectively detected by the PMT, accounting for its limited angular acceptance.

What I have decided to do is divide reflectance computation by this fraction here, to account for the fact that the reflected PMT sees less light.

# 25/08/2025
In our setup, we compute the fraction of reflected photons hitting the PMT as:
$$
f_{\text{geom}} = \frac{\int_{\theta_{\min}}^{\theta_{\max}} F(\theta)  d\theta}{\int_{-90^\circ}^{90^\circ} F(\theta)  d\theta}$$
where \(F(\theta)\) is the **angular reflection function** derived from the PTFE paper.

However, the PMT is **circular**, not square. This introduces a non-uniform acceptance across the azimuthal angle: photons at larger polar angles hit the circular hole with different probabilities. To account for this, we define an effective acceptance correction \(C(x)\) as:
$$
C(x) = \sqrt{1 - \left(\frac{x-a}{R}\right)^2}$$
where:  
- \(R = 25^\circ\) represents **half of the PMT angular acceptance**,  
- \(a\) is the incidence angle of the PMT,  
- \(x\) is the azimuthal angle of the photon.

## Computation of \(C(x)\)
The factor \(C(x)\) is derived from the **geometry of a circle**: it is the ratio between the length of the chord at a given distance from the center and the full diameter of the circle. For a circular PMT:
$$
\text{chord length} = 2 \sqrt{R^2 - (x-a)^2} \quad \Rightarrow \quad C(x) = \frac{\text{chord length}}{2R} = \sqrt{1 - \left(\frac{x-a}{R}\right)^2}$$
This captures how the **effective acceptance decreases** for photons hitting away from the center of the PMT.
The corrected transmission function is then:
$$
T(x) = F(x) \cdot C(x)$$
so that the geometrically weighted fraction becomes:
$$
f_{\text{geom}} = \frac{\int_{x_{min}}^{x_{max}} F(x) \cdot C(x)  dx}{\int_{-90^\circ}^{90^\circ} F(x)  dx}$$

This approach effectively accounts for the **circular geometry of the PMT**, providing a realistic estimate of the photons that contribute to the PMT signal.

# 27/08/2025
The "transmission coefficient" or effective attenuation length (λ) you see in papers for PTFE is not a fundamental constant. It varies significantly because it's an effective property that depends heavily on both the specific physical structure of the PTFE sample and its surrounding environment.

Here’s why the values are different, especially when used with xenon:

## 1. It's Not Really the Beer-Lambert Law
First, let's clarify the physics. While researchers might use the term for convenience, the Beer-Lambert law  
\( I = I_0 e^{-\alpha z} \) technically only applies to absorption in a transparent medium.  

In PTFE, the light intensity drops primarily because of intense scattering.

So, the coefficient they are measuring is an **effective attenuation length** \( L_{\text{att}} \).  
It describes the combined effect of some absorption \( \mu_a \) and a lot of scattering \( \mu_s \).  
It's a practical parameter that tells you how far light generally penetrates, but it's not a true absorption coefficient.

## 2. The Surrounding Medium is Critical (Xenon vs. Air)
This is the most important factor based on your question. PTFE is a porous, sintered material filled with microscopic voids. The scattering happens at the interface between the PTFE polymer and whatever is filling these voids.

The amount of scattering depends on the mismatch in the index of refraction (n) at this interface.

- **In Air**:  
  The voids are filled with air (\( n_{\text{air}} \approx 1.0 \)).  
  The refractive index of PTFE is \( n_{\text{PTFE}} \approx 1.35 \).  
  The mismatch (1.35 vs. 1.0) is large, causing very strong scattering.  
  → High reflectivity, shorter attenuation length.

- **In Liquid Xenon (LXe)**:  
  The voids become filled with liquid xenon (\( n_{\text{LXe}} \approx 1.6 \)).  
  The mismatch is now between PTFE (1.35) and LXe (1.6).  
  → Reduced scattering, longer attenuation length.

- **In Gaseous Xenon (GXe)**:  
  The gas will have an index higher than air, again changing the scattering conditions.

Essentially, when you immerse PTFE in a liquid or dense gas, you change what's inside its pores, which fundamentally alters how it scatters light.  
A medium with a refractive index closer to that of PTFE will reduce scattering and make the material appear more transparent, increasing the attenuation length.

## 3. Not All PTFE is Created Equal
Even if tested in the same medium, the PTFE itself can vary significantly:

- **Density and Porosity**:  
  The primary factor. PTFE is made by sintering (compressing and heating) a powder.  
  Higher sintering pressure → denser material → fewer/smaller voids.  
  Denser PTFE has less internal surface area for scattering → generally increases attenuation length.  
  (Papers often report the density of their PTFE, e.g. g/cm³, for this reason.)

- **Manufacturing Process**:  
  The brand (e.g., *Teflon®*, *Spectralon®*), preparation (sintered, extruded), and surface finish (machined smooth, left rough) all change the microstructure and, therefore, the scattering properties.

- **Purity**:  
  Industrial-grade PTFE may contain contaminants that absorb light, decreasing the measured transmission.

## In Summary
You're seeing different values because researchers are measuring **different physical systems**.  
The "PTFE" is just one component. The complete system includes:

- the material's specific density,  
- its purity,  
- and the liquid or gas it's immersed in.  

Each combination will yield a different, but valid, effective attenuation length.

# 28/08/2025
# Reflectance Correction and Error Propagation

I've added error propagation for the reflectance correction. The sources of uncertainty are:

1. **Poissonian errors** in the number of photoelectrons \(N_\text{PE}\), computed as \(\sqrt{N_\text{PE}}\).  
2. **Geometrical errors** associated with the parameter \(a = 25^\circ\) in the `C(x)` function, representing half of the PMT angular acceptance.  
3. **Angular uncertainty** in the angle of incidence/reflection, \(\theta_\text{inc} = \theta_\text{refl}\).  
4. **Parameter uncertainties** in the fit function \(F(x)\): \(R_1\), \(R_2\), and \(\sigma\).

## Definition of the Correction Factor

The transmittance correction factor is defined as:
$$
\text{corrT} = \frac{1}{D} \int_{\theta_1}^{\theta_2} T(x)  dx, 
\quad
T(x) = F(x) \cdot C(x)
$$
where:
$$
D = \int_{-90}^{90} F(x)  dx, \quad
C(x) = \sqrt{1 - \left(\frac{x - a}{R}\right)^2}, \quad
\theta_\text{refl} = \theta_\text{inc}, \quad a = 25.
$$

## Linear Error Propagation

We apply **linear propagation** (first-order approximation) for the uncertainties in \(a\), \(\theta_\text{refl}\), \(R_1\), \(R_2\), and \(\sigma\):
$$
\delta(\text{corrT}) \approx \frac{1}{D} \left[
\delta a \int_{\theta_1}^{\theta_2} \frac{\partial T(x)}{\partial a} dx
+ \delta \theta_\text{refl} \int_{\theta_1}^{\theta_2} \frac{\partial T(x)}{\partial \theta_\text{refl}} dx
+ \delta R_1 \int_{\theta_1}^{\theta_2} \frac{\partial T(x)}{\partial R_1} dx
+ \delta R_2 \int_{\theta_1}^{\theta_2} \frac{\partial T(x)}{\partial R_2} dx
+ \delta \sigma \int_{\theta_1}^{\theta_2} \frac{\partial T(x)}{\partial \sigma} dx
\right]
$$

### Partial derivatives

1. **Derivative w.r.t \(a\) (PMT acceptance):**
$$
\frac{\partial T(x)}{\partial a} = F(x) \cdot \frac{\partial C(x)}{\partial a}, 
\quad
\frac{\partial C(x)}{\partial a} = \frac{(x - a)^2}{a^3 \sqrt{1 - ((x - a)/a)^2}}
$$

2. **Derivative w.r.t \(\theta_\text{refl}\) (angle of incidence):**
$$
\frac{\partial T(x)}{\partial \theta_\text{refl}} = F(x) \cdot \frac{\partial C(x)}{\partial \theta_\text{refl}}, 
\quad
\frac{\partial C(x)}{\partial \theta_\text{refl}} = \frac{(x - \theta_\text{refl})}{a^2 \sqrt{1 - ((x - \theta_\text{refl})/a)^2}}
$$

3. **Derivative w.r.t \(R_1\), \(R_2\), and \(\sigma\) (fit parameters):**

- **Derivative with respect to \(R_1\):**
$$
\frac{\partial T(x)}{\partial R_1} = C(x) \cdot \frac{\partial F(x)}{\partial R_1}, \quad
\frac{\partial F(x)}{\partial R_1} = \frac{1}{\text{Norm}} \exp\Bigg[-\frac{\alpha(x, \theta)^2}{2 \sigma^2}\Bigg]
$$

- **Derivative with respect to \(R_2\):**
$$
\frac{\partial T(x)}{\partial R_2} = C(x) \cdot \frac{\partial F(x)}{\partial R_2}, \quad
\frac{\partial F(x)}{\partial R_2} = \frac{1}{\text{Norm}} \cos\theta
$$

- **Derivative with respect to \(\sigma\):**
$$
\frac{\partial T(x)}{\partial \sigma} = C(x) \cdot \frac{\partial F(x)}{\partial \sigma}, \quad
\frac{\partial F(x)}{\partial \sigma} = \frac{R_1 \exp[-\alpha(x,\theta)^2 / (2\sigma^2)] \, \alpha(x,\theta)^2}{\sigma^3 \, \text{Norm}}
$$

where:  

- \(C(x)\) is the geometrical acceptance function.  
- \(\alpha(x, \theta) = \alphaFunc(x, \theta)\) is the argument inside the exponential of the fit function.  
- \(\text{Norm}\) is the normalization factor of the fit.  
- \(\theta\) is the angle of incidence/reflection.

# Computation of the Average Pulse Rate

When computing the **rate of pulses** in a given time region, it is important to normalize correctly by the total **observation time**.

## 1. Simple rate per waveform

If you just sum the pulse areas in a region of a single waveform, the rate would be:
$$
R_\text{simple} = \frac{\text{Total pulse area in region}}{\Delta t_\text{region}}
$$

where:

- \(\Delta t_\text{region}\) is the duration of the time window (e.g., pre-trigger, trigger, post-trigger).  

**Problem:** This ignores how many waveforms you analyzed and overestimates the true average rate.

## 2. Corrected average rate

The **correct rate** should account for all the waveforms analyzed. If you have:

- \(\Delta t_\text{region}\) = duration of the time region,  
- \(N_\text{waveforms}\) = number of waveforms analyzed (including those without pulses in the region),  
- \(\text{Area}_i\) = total pulse area in the region for the \(i\)-th waveform,  

then the **average rate** is:
$$
R_\text{avg} = \frac{\sum_{i=1}^{N_\text{waveforms}} \text{Area}_i}{\Delta t_\text{region} \cdot N_\text{waveforms}}
$$

### 3. Explanation

1. **Numerator:** sum of pulse areas only in the chosen region (e.g., trigger region) across all waveforms.  
2. **Denominator:** total “lifetime” or total observation time in that region:
$$
t_\text{total} = \Delta t_\text{region} \cdot N_\text{waveforms}
$$

3. **Why include all waveforms:**  
   - Waveforms without pulses still count as “observed time” during which no signal occurred.  
   - Excluding them would **overestimate the rate**, because you would ignore the silent periods.

### 4. Summary

- **Step 1:** For each waveform, sum the areas of pulses in the region.  
- **Step 2:** Sum these areas over all waveforms.  
- **Step 3:** Divide by the total observation time (\(\Delta t_\text{region} \times N_\text{waveforms}\)) to get the average rate:
$$
\boxed{R_\text{avg} = \frac{\sum_{i=1}^{N_\text{waveforms}} \text{Area}_i}{\Delta t_\text{region} \cdot N_\text{waveforms}}}
$$

This gives the **true average pulse rate per unit time**, including both waveforms with and without pulses.
