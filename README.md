# Build-Your-Own-Radar-System-Project
This project involved building two practical radar systems - a short-range circuit-based radar operating at 2.4 GHz, and a software-defined radar implementation at 5.8 GHz using GNU Radio.

The goal was to gain hands-on experience with radar fundamentals, more specifically:

* Circuit implementation for short range radar at 2.4 GHz
* Range measurement using an FMCW radar
* Velocity measurement employing a CW radar
* Landscape imaging exploiting SAR technique
* Implementing a radar system using a software defined radio at 5.8GHz.
  
# Circuit implementation
The circuit radar utilized an voltage controlled oscillator driven by a modulator chip to generate continuous wave and FMCW signals. A custom video amplifier circuit was built to condition receive signals.

# Software Defined Radar
The SDR was implemented using a USRP B205mini and GNU Radio. Block diagrams were developed for continuous wave and low-IF CW modes for both range and velocity measurement.

# Signal Processing
All the matlab codes are implemented for:

* CW velocity detection through FFT algorithms
* FMCW range measurement using the IFFT
* Range Migration Algorithm for SAR signal processing
* QM signal processing for velocity direction suing the SDR.

# Results
Both systems were tested and able to capture:

* Single/multiple target range and velocity
* Landscape imaging
* Physiological signals like breathing and heart rate
