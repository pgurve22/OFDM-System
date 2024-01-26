# OFDM-System

## Single-Input-Single-Output OFDM System


![ofdm_blk](ofdm_blk.png)


- Problem Statement : Simulate an OFDM wireless system in MATLAB with 𝑁 = 64 subcarriers for a channel with 𝐿 = 3 i.i.d. Rayleigh fading unit gain channel taps. Generate the BER curves vs dB SNR for BPSK symbols loaded over all the subcarriers and also superimpose the plots obtained via the corresponding analytical expression derived in class lectures. Choose the SNR range so as to obtain BER values up to 10−4.
  
- Simulation :
  

![SISO_OFDM_BER](SISO_OFDM_BER.png)


## Single-Input-Multi-Output OFDM System

- Problem Statement : Extend the above simulation to a SIMO OFDM system with 𝑅 = 2 receive antennas and similarly plot the BER curves obtained via both simulation and analysis. Submit the code and relevant plots for both the problems above.
 
- Simulation :
  
![SIMO_OFDM_BER](SIMO_OFDM_BER.png)


## Multi-Input-Multi-Output OFDM System

### Transmitter 


![mimo_ofdm_tx](mimo_ofdm_tx.png)


### Receiver


![mimo_ofdm_rx](mimo_ofdm_rx.png)


- Problem Statement : Simulate a MIMO OFDM wireless system in MATLAB with 𝑁 = 64 subcarriers for a channel with 𝐿 = 3 i.i.d. Rayleigh fading unit gain channel taps, 𝑡 = 2 transmit and 𝑟 = 2 receive antennas. Generate the BER curves vs dB SNR for BPSK symbols loaded over all the subcarriers and all transmit antennas. Superimpose the plots obtained via the corresponding analytical expression derived in class lectures. Choose the SNR range so as to obtain BER values up to 10−4
 
- Simulation :
  
![MIMO_OFDM_BER](MIMO_OFDM_BER.png)


## MU-MIMO-OFDM Downlink System

- Problem Statement : Simulate a MU-MIMO-OFDM wireless system in MATLAB with M = 32 subcarriers for a channel with 𝐿 = 5 i.i.d. Rayleigh fading unit gain channel taps and Pedestrian-A channel with gain [0 −9.7 −19.2 −22.8], N_t = 64 BS Antennas and K = 4 Single Antenna Users. Generate the both Uplink and Downlink sum rate for QAM symbols loaded over all the subcarriers and compared it for MRC,ZF,MMSE.
 
- Simulation : Rayleigh Fading Channel 


![MU_MIMO_OFDM_SE_DL](Downlink_Sum_Rate_Ray.png)


- Simulation : Pedestrian-A Channel 


![MU_MIMO_OFDM_SE_DL](Downlink_Sum_Rate_PedA.png)


## MU-MIMO-OFDM Uplink System

- Simulation : Rayleigh Fading Channel 


![MU_MIMO_OFDM_SE_UL](Uplink_Sum_Rate_Ray.png)


- Simulation : Pedestrian-A Channel 


![MU_MIMO_OFDM_SE_UL](Uplink_Sum_Rate_PedA.png)
