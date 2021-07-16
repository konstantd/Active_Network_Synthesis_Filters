# Active_Network_Synthesis_Filters



This repository contains the implementation of four different filters:



- High-Pass Butterworth
- Low-Pass Chebyshev 
- Band-Elimination Inverse Chebyshev 
- Band-Pass Chebyshev 



Originally, the analytical solutions for each filter was implemented in MATLAB. Also, MATLAB was used for the calculation of the parameters of the electronic components. After that, the schematic was designed in MULTISIM and the application was compared with the analytical solution in order to assure that the respective filter works properly. The scripts, multisim files as well as the pdf files with detailed description of the whole steps that were folllowed can be found in the respective folders for each filter. 

Below, we are going to present the requirements, the final Transfer Function with all Units that were used, the schematic , a Transient Analysis as well as a Fourier Analysis in Multisim of the input and output signal, in order to prove the efficiency of the designed filter.



### High-Pass Butterworth



We would like to implement a Band Elimination Inverse Chebyshev filter with the below requirements of frequencies and attenuation levels:

- fp = 5 kHz 
- fs = 1.92 kHz 
- amax = 0.5278 dB 
- amin = 29.3333 dB
- 10 dB Gain



In the below graph we can see the Transfer function of the filter. It satisfies all the above requirements. 

<p allign = "center">
     <img src="/HP_Butterworth/photos/HP_1.png"width = "80%">
</p>

Proceeding with the circuit, after using the script to calculate all the required components, below we can see the final schematic.

<p allign = "center">
     <img src="/HP_Butterworth/photos/HP_2.png"width = "80%">
</p>

Connecting a Bode-Plotter, we take exactly the same above Transfer function graph in Multisim. We are confident that the filter works properly.

<p allign = "center">
     <img src="/HP_Butterworth/photos/HP_3.png"width = "80%">
</p>

We are going to test the filter, giving as an input a sum of different cosines using different AC voltage sources in series and an Oscilloscope. The output of the Transient Analysis with the Input (green) and Output (red) signals can be seen below. 

<p allign = "center">
     <img src="/HP_Butterworth/photos/HP_4.png"width = "80%">
</p>

Proceeding with the Fourier Analysis, below we can see the spectrum of the input signal:

<p allign = "center">
     <img src="/HP_Butterworth/photos/HP_5.png"width = "80%">
</p>

Spectrum of the output signal:

<p allign = "center">
     <img src="/HP_Butterworth/photos/HP_6.png"width = "80%">
</p>

We observe that low frequencies (below fs = 1.9231 kHz) are "cut" in the output. At the input we have 5 spikes while the output lacks the frequencies for 384.6 Hz and for 1346 Hz that are indeed below fs, which makes perfect sense since our circuit is an High pass filter. At the same time, the correctness of the gain regulation is observed. The magnitude of the spikes at the output is greater (by approximately 3.1 times) than the spikes at the input, since we have a gain of 10dB. We conclude that the filter works properly, as all its requirements are met.



### Low-Pass Chebyshev 



We would like to implement a Band Elimination Inverse Chebyshev filter with the below requirements of frequencies and attenuation levels:

- fp = 6 kHz
-  fs = 10.320 kHz 
- amax = 0.25 dB 
- amin = 26.25 dB 
- 0 dB

In the below graph we can see the Transfer function of the filter. It satisfies all the above requirements. 

<p allign = "center">
     <img src="/LP_Chebyshev/photos/LP_1.png"width = "80%">
</p>

Proceeding with the circuit, after using the script to calculate all the required components, below we can see the final schematic.

<p allign = "center">
     <img src="/LP_Chebyshev/photos/LP_2.png"width = "80%">
</p>

Connecting a Bode-Plotter, we take exactly the same above Transfer function graph in Multisim. We are confident that the filter works properly.

<p allign = "center">
     <img src="/LP_Chebyshev/photos/LP_3.png"width = "80%">
</p>

We are going to test the filter, giving as an input a rectangular pulse. The output of the Transient Analysis with the Input (green) and Output (red) signals can be seen below. 

<p allign = "center">
     <img src="/LP_Chebyshev/photos/LP_4.png"width = "80%">
</p>

Proceeding with the Fourier Analysis, below we can see the spectrum of the input signal:

<p allign = "center">
     <img src="/LP_Chebyshev/photos/LP_5.png"width = "80%">
</p>

Spectrum of the output signal:

<p allign = "center">
     <img src="/LP_Chebyshev/photos/LP_6.png"width = "80%">
</p>

We observe that high frequencies above fs, ie above 10.32 kHz are "cut" in the output, and hence these frequencies that create the rectangular pulse at the input signal have zero magnitude in the output range, which makes sense since our circuit is a low-pass filter. At the same time, the correctness of the gain regulation is observed. The magnitude of the spikes at the output for low frequencies is almost the same as the input since we have a 0 dB gain. So it is concluded that the filter works properly as all requirements are satisfied.

 

### Band Elimination Inverse Chebyshev

We would like to implement a Band Elimination Inverse Chebyshev filter with the below requirements of frequencies and attenuation levels:

- f0 = 2.4 kHz		

- f1 = 1.75 kHz

-  f2 = 3.2914 kHz
-  f3 = 2.0751 kHz
-  f4 = 2.7758 kHz 
- αmin = 32.4444 dB
- αmax = 0.5556 dB
- 10 dB gain



In the below graph we can see the Transfer function of the filter. It satisfies all the above requirements. 

<p allign = "center">
     <img src="/BE_Inverse_Chebyshev/photos/BE_1.png"width = "80%">
</p>

Proceeding with the circuit, after using the script to calculate all the required components, below we can see the final schematic.

<p allign = "center">
     <img src="/BE_Inverse_Chebyshev/photos/BE_1.png"width = "80%">
</p>

Connecting a Bode-Plotter, we take exactly the same above Transfer function graph in Multisim. We are confident that the filter works properly.

<p allign = "center">
     <img src="/BE_Inverse_Chebyshev/photos/BE_6.png"width = "80%">
</p>

We are going to test the filter, giving as an input a sum of different cosines using different AC voltage sources in series and an Oscilloscope. The output of the Transient Analysis with the Input (green) and Output (red) signals can be seen below. 

<p allign = "center">
     <img src="/BE_Inverse_Chebyshev/photos/BE_3.png"width = "80%">
</p>

Proceeding with the Fourier Analysis, below we can see the spectrum of the input signal:

<p allign = "center">
     <img src="/BE_Inverse_Chebyshev/photos/BE_4.png"width = "80%">
</p>

Spectrum of the output signal:

<p allign = "center">
     <img src="/BE_Inverse_Chebyshev/photos/BE_5.png"width = "80%">
</p>



At the input of the filter there are five spikes while at the output of one of them,  f = 2.2376 kHz (which is indeed the only one frequency between f3 = 2.0751 kHz and f4 = 2775.76 kHz) is eliminated, something that is expected since our circuit is a band elimination filter. At the same time, the correctness of the gain is observed. The width of the spikes at the output are greater (by about 3 times) than the spikes at the input, since we have a gain of 10dB. This concludes that the filter is working properly, as all requirements are met.



### Band-Pass Chebyshev 

We would like to implement a Band Elimination Inverse Chebyshev filter with the below requirements of frequencies and attenuation levels:

- f0 = 900 Hz
- f1 = 850 Hz
-  f2 = 952.941 Hz
- f3 = 793.8602 Hz
- f4 = 1.0203 kHz 
- αmin= 29.0556 dB
- αmax = 0.7222 dB		
- 0 dB Gain



In the below graph we can see the Transfer function of the filter. It satisfies all the above requirements. 

<p allign = "center">
     <img src="/BP_Chebyshev/photos/BP_1.png"width = "80%">
</p>

Proceeding with the circuit, after using the script to calculate all the required components, below we can see the final schematic.

<p allign = "center">
     <img src="/BP_Chebyshev/photos/BP_2.png"width = "80%">
</p>

Connecting a Bode-Plotter, we take exactly the same above Transfer function graph in Multisim. We are confident that the filter works properly.

<p allign = "center">
     <img src="/BP_Chebyshev/photos/BP_3.png"width = "80%">
</p>

We are going to test the filter, giving as an input a sum of different cosines using different AC voltage sources in series and an Oscilloscope. The output of the Transient Analysis with the Input (green) and Output (red) signals can be seen below. 

<p allign = "center">
     <img src="/BP_Chebyshev/photos/BP_4.png"width = "80%">
</p>

Proceeding with the Fourier Analysis, below we can see the spectrum of the input signal:

<p allign = "center">
     <img src="/BP_Chebyshev/photos/BP_5.png"width = "80%">
</p>

Spectrum of the output signal:

<p allign = "center">
     <img src="/BP_Chebyshev/photos/BP_6.png"width = "80%">
</p>

Initially, we noticet frequencies 𝑓2 = 1483.33 𝐻𝑧, 𝑓3 = 317.54 𝐻𝑧, 𝑓4 = 2.55 𝑘𝐻𝑧, 𝑓5 = 3.06 𝑘𝐻𝑧, which are in the cut-off zone are "cut", ie they are eliminated in the output . At the output we have only one spike at the frequency of 𝑓1 = 875 𝐻𝑧 which is actually between 850 and 952.94 Hz, which it absolutely logical since our circuit is a bandpass filter. Also, the correctness of the gain regulation is also observed. The magnitude of the spikes remaining at the output is approximately equal to its amplitude at the input, since we have a gain of 0dB. This concludes that the filter works correctly, as all requirements are met.

