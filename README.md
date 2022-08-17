# DAGAF

Lin, Yue-Der, Yong Kok Tan, and Baofeng Tian. "A novel approach for decomposition of biomedical signals in different applications based on data-adaptive Gaussian average filtering." Biomedical Signal Processing and Control 71 (2022): 103104.

DOI = https://doi.org/10.1016/j.bspc.2021.103104

## Abstract
The analysis of biomedical signals plays a crucial role in modern medicine and physiological research. Due to
exogenous or endogenous interferences and complex interaction among physiological systems, biomedical signals
usually possess nonstationary characteristics. To extract the features buried in such kind of signals, the signal
decomposition algorithm that is data-adaptive and stable is highly demanded. This paper introduces a novel
decomposition algorithm termed as data-adaptive Gaussian average filtering (DAGAF) that is new and potential
in biomedical applications. Six biomedical scenarios, including finger photoplethysmography (PPG) signal with
obscure respiratory-induced intensity variation (RIIV) component, wrist PPG signal with apparent RIIV
component, seismocardiography (SCG) signal with implicit respiration component, electrocardiogram (ECG)
with baseline wander (BW), ECG with power-line interference (PLI), and R-R intervals (RRI) sequence with
implicit low-frequency trend wave, are adopted as examples for computer experiments. The results of computer
experiments verify that the DAGAF algorithm can satisfy the specified requirement in different biomedical
scenarios. DAGAF algorithm possesses the advantages of mathematical formulation and computational efficiency.
It can be an alternative choice besides the commonly used empirical mode decomposition (EMD) and
ensemble empirical mode decomposition (EEMD).

## DAGAF in Python

### Requirements

1. Python >= 2.7
- Install Python using the Anaconda is reommand.
2. numpy
- Install using pip (``pip install numpy``)
3. scipy
- install using pip (``pip install scipy``)
4. matplotlib
- install using pip (``pip install matplotlib``)

### Usage
An simple example is as follows.
```
import numpy as np
import matplotlib.pyplot as plt
from dagaf import dagaf
s = np.random.random(1000)

# Plot the data.
plt.figure()
plt.plot(s)
plt.show()
imf, res = dagaf(s, 4, 1.6, 'd', 'd', 20, 0.001, True)
```
The results would show three figures:

<p align="left">
  1. The computation steps of DAGAF <br />
  <img src="https://github.com/t22302856/DAGAF/blob/main/Fig_Signal_Ave.jpg" width="450" height="300" alt="Signal_Ave"/> <br />
  2. The IMF signals <br />
  <img src="https://github.com/t22302856/DAGAF/blob/main/Fig_imf.jpg" wwidth="450" height="300" alt="imf"/> <br />
  3. The residual component after DAGAF <br />
  <img src="https://github.com/t22302856/DAGAF/blob/main/Fig_Signal_Res.jpg" width="450" height="300" alt="Signal_Res"/>  
</p>

## DAGAF in Matlab

### Running environment

- MATLAB 2020b

### Usage
A simple example for Matlab.

```
s = rand(1000,1);
[imf, res] = dagaf(s, 4, 1.6, 'd', 'd', 20, 0.001);
```

## Citation
If you would like to use this package and function, please cite this paper:

```
@article{LIN2022103104,
title = {A novel approach for decomposition of biomedical signals in different applications based on data-adaptive Gaussian average filtering},
journal = {Biomedical Signal Processing and Control},
volume = {71},
pages = {103104},
year = {2022},
issn = {1746-8094},
doi = {https://doi.org/10.1016/j.bspc.2021.103104},
url = {https://www.sciencedirect.com/science/article/pii/S1746809421007011},
author = {Yue-Der Lin and Yong Kok Tan and Baofeng Tian},
keywords = {Baseline wander (BW), Data-adaptive Gaussian average filtering (DAGAF), Power-line interference (PLI), Respiratory-induced intensity variation (RIIV), Signal decomposition},
abstract = {The analysis of biomedical signals plays a crucial role in modern medicine and physiological research. Due to exogenous or endogenous interferences and complex interaction among physiological systems, biomedical signals usually possess nonstationary characteristics. To extract the features buried in such kind of signals, the signal decomposition algorithm that is data-adaptive and stable is highly demanded. This paper introduces a novel decomposition algorithm termed as data-adaptive Gaussian average filtering (DAGAF) that is new and potential in biomedical applications. Six biomedical scenarios, including finger photoplethysmography (PPG) signal with obscure respiratory-induced intensity variation (RIIV) component, wrist PPG signal with apparent RIIV component, seismocardiography (SCG) signal with implicit respiration component, electrocardiogram (ECG) with baseline wander (BW), ECG with power-line interference (PLI), and R-R intervals (RRI) sequence with implicit low-frequency trend wave, are adopted as examples for computer experiments. The results of computer experiments verify that the DAGAF algorithm can satisfy the specified requirement in different biomedical scenarios. DAGAF algorithm possesses the advantages of mathematical formulation and computational efficiency. It can be an alternative choice besides the commonly used empirical mode decomposition (EMD) and ensemble empirical mode decomposition (EEMD).}
}
```

## Contact
Feel free to contact Kai-Chun Liu (t22302856@gmail.com) or Yu-Der Lin (ydlin@fcu.edu.tw). Welcom any quesntion or discussion. 
