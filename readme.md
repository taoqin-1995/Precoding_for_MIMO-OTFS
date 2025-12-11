# Precoding_for_MIMO-OTFS：Capacity Analysis of Single-User MIMO-OTFS Systems

This repo is the MATLAB simulation code for the paper "Capacity Analysis of Single-User MIMO-OTFS Systems".

Q. Tao and W. Yuan, "Capacity Analysis of Single-User MIMO-OTFS Systems," in IEEE Transactions on Vehicular Technology, doi: 10.1109/TVT.2025.3588218.

https://ieeexplore.ieee.org/document/11078311 

## **Paper Overview**

This paper analyzes the capacity of single-user MIMO-OTFS systems and provides corresponding transceiver designs, including precoding and equalization schemes.

In fact, the analytical framework established in this work is readily extensible to other emerging waveform systems, such as MIMO-AFDM... 

 

## **Repository Structure**

##### 【Capacity_Asymptotic_capacity.m】

Evaluates exact and asymptotic capacity performance and reproduces Fig. 2 of the paper.

<div align=center>
<img src="./figures/exact and asymptotic capacity performance.png" width = "550" height = "450" alt="">
</div>

##### 【Capacity_beamforming.m】

Simulates capacity under different beamforming schemes and reproduces Fig. 3 of the paper.

<div align=center>
<img src="./figures/capacity under different beamforming schemes.png" width = "550" height = "450" alt="">
</div>

## Citation

If this code helps you, please consider citing our works:

```
@article{tao2025capacity,
 author={Tao, Q. and Yuan, W.},
 journal={IEEE Transactions on Vehicular Technology},
 title={Capacity Analysis of Single-User MIMO-OTFS Systems},
 year={2025},
 doi={10.1109/TVT.2025.3588218}
}
```



 