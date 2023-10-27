# Build-Your-Own-Radar-System-Project
This project involved building two practical radar systems - a short-range circuit-based radar operating at 2.4 GHz, and a software-defined radar implementation at 5.8 GHz using GNU Radio.

The goal was to gain hands-on experience with radar fundamentals, more specifically:

* Circuit implementation for short range radar at 2.4 GHz
* Range measurement using an FMCW radar
* Velocity measurement employing a CW radar
* Landscape imaging exploiting SAR technique
* Implementing a radar system using a software defined radio at 5.8GHz.
  
# Circuit implementation for short range randar at 2.4Ghz
The circuit radar utilized an voltage controlled oscillator driven by a modulator chip to generate continuous wave and FMCW signals. A custom video amplifier circuit was built to condition receive signals.
![Circuit overall](https://github.com/GiandomenicoPanettieri/Build-Your-Own-Radar-System-Project/assets/57620194/03967a01-5e8d-4ce9-b239-fee25021f8dd)

# Signal Processing
All the matlab codes are implemented for:

* CW velocity detection through FFT algorithms
* FMCW range measurement using the IFFT
* Range Migration Algorithm for SAR signal processing
* QM signal processing for velocity direction suing the SDR.

![immagine_2023-10-27_142351985](https://github.com/GiandomenicoPanettieri/Build-Your-Own-Radar-System-Project/assets/57620194/6f2e410c-7060-4c8b-8c85-1a124a53ba95)

# Software Defined Radar system project at 5.8GHz
The SDR was implemented using a USRP B205mini and GNU Radio. Block diagrams were developed for continuous wave and low-IF CW modes for both range and velocity measurement.

# Results
Both systems were tested and able to capture:

* Single/multiple target range and velocity
* Landscape imaging
* Physiological signals like breathing and heart rate
![BreathingMarko](https://github.com/GiandomenicoPanettieri/Build-Your-Own-Radar-System-Project/assets/57620194/4d9abf06-1376-441d-b69f-89d60db97f06)

