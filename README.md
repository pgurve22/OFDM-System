# OFDM-System

## Single-Input-Single-Output OFDM System


Total number of sub-carriers is $M$. Hence total number of information symbols is also $M$.


$$X(0) , X(1) , X(2) ,\space.\space .\space .\space.\space ,X(M-1)$$


So basically $u$-th IFFT coefficient of the information symbols $X(0) , X(1) , X(2) ,\space.\space .\space .\space.\space, X(M-1)$ is given as-


```math
    x(u) = \underbrace{\sum_{i} X(i) e^{2\pi\frac{iu}{M}}}_{IFFT}
```

Similarly $x(0),x(1),\space.\space .\space .\space.\space,x(M)$ are the IFFT coefficient at the $0,1,\space.\space .\space .\space.\space,M$ time points respectively. Our serial signal is given as-


```math
    \hat{x} = [x(0),x(1),x(2),\space.\space .\space .\space.\space,x(M)]
```


- Cyclic Prefix in OFDM -

Consider a frequency-selective channel modelled with channel taps $h(0),h(1),\space.\space .\space .\space.\space,h(L_h-1)$. Thus, the received symbol $y$ at a given time instant $n$ can be
expressed as  


```math
    y(n)=h(0)x(n)+\underbrace{h(1)x(n-1)+h(2)x(n-2)+\space.\space .\space .\space.\space+h(L_h-1)x(n-L_h+1)}_{ISI \space Component}
```


Consider now two OFDM symbols as follows. Let $x (0)$ , $x (1)$ , . . . . ., $x (N − 1)$ denote the IFFT samples of the modulated symbols $X(0)$ , $X (1)$ , . . . . . , $X (N − 1)$, while $\tilde{x}(0)$ , $\tilde{x}(1)$ , . . . . . , $\tilde{x}(N − 1)$ denote the IFFT samples of the previous modulated symbol block $\tilde{X}(0)$ , $\tilde{X}(1)$ , . . . . . , $\tilde{X}(N − 1)$. Thus, the samples corresponding to these two blocks of OFDM symbols are transmitted sequentially as


```math
    \underbrace{\tilde{x}(0) , \tilde{x}(1) , \space.\space .\space . \space, \tilde{x}(N − 1}_{Previous \space block} \space, \space \underbrace{x (0) , x (1) , \space.\space .\space .\space , x (N − 1)}_{Current \space block}
```


```math
    y(0)=h(0)x(0)+\underbrace{h(1)\tilde{x}(N-1)\space+\space .\space .\space .\space.\space\space+\space h(L_h-1)\tilde{x}(N-L_h+1)}_{ISI \space from \space the \space previous \space OFDM \space symbol}
```


```math
    y(1)=h(0)x(1)+h(1)x(0)+\underbrace{h(2)\tilde{x}(N-1)\space+\space .\space .\space .\space.\space\space+\space h(L_h-1)\tilde{x}(N-L_h+2)}_{ISI \space from \space the \space previous \space OFDM \space symbol}
```


Let us now consider a modified transmission scheme as follows. To each transmitted OFDM sample stream, we pad the last $L_c$ symbols to make the transmitted stream as follows.


```math
    \underbrace{\tilde{x}(0) , \tilde{x}(1) , \space.\space .\space . \space, \tilde{x}(N − 1}_{Previous \space block} \space, \space \underbrace{x (N-L_c) , x (N-L_c+1) , \space.\space .\space .\space , x (N − 1)}_{Cyclic \space prefix} \space, \space \underbrace{x (0) , x (1) , \space.\space .\space .\space , x (N − 1)}_{Current \space block}
```


Consider now the received symbol corresponding to $x (0)$. This is given as


```math
    y(0)=h(0)x(0)+\underbrace{h(1)x(N-1)\space+\space .\space .\space .\space.\space\space+\space h(L_h-1)x(N-L_h+1)}_{ISI \space from \space the \space same \space OFDM \space symbol}
```


The inter-symbol interference can be seen to now be from $x (N − 1) , x (N − 2) , . . . , x (N − L + 1)$, if $L_c ≥ L − 1$. Thus, with the cyclic prefix of appropriate length, i.e., $L_c ≥ L − 1$, inter-OFDM symbol interference can be avoided and inter-symbol interference is restricted to samples from the same OFDM symbol. Therefore, the samples $y (0) , y (1) , . . . , y (N − 1)$ are given as



```math
    y (0) = h (0) x (0) \space+ \space h (1) x (N − 1) \space+ \space .\space . \space . \space+ \space h (L − 1) x (N − L + 1)   
```


```math
    y (1) = h (0) x (1) \space + \space h (1) x (0) \space+ \space \space .\space . \space . \space+ \space h (L − 1) x (N − L + 2)
```


```math
                            . \space . \space . \space . \space .
```


```math
    y (N − 1) = h (0) x (N − 1) \space + \space h (1) x (N − 2) \space + \space .\space . \space . \space + \space h (L − 1) x (N − L)
```


It can now be clearly seen that the output $y (n)$ is a circular convolution between the channel filter $h (n)$ and the input $x (n)$.



```math
    [y (0) \space,\space y (1)  \space,\space . \space . \space .  \space,\space y (N − 1)] = [h (0)  \space,\space h (1)  \space,\space . \space . \space .  \space,\space h (L − 1)  \space,\space 0  \space,\space . \space . \space .  \space,\space 0] \ast_N \space [x (0)  \space,\space
x (1)  \space,\space . \space . \space .  \space,\space x (N − 1)]
```


where $\ast_N$ denotes circular convolution of modulo $N$. Therefore, the output $y$ can be written as-


```math
    y = h \ast_N x
```


Therefore, taking the FFT of $y (n)$ at the output, we have


```math
    Y (k) = H (k) \space X (k) \space \space \space , \space \space \space 0 ≤ k ≤ N − 1
```






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
