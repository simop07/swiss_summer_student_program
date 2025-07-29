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