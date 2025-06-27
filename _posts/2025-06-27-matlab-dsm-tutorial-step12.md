---
title: How to design NTF using MATLAB
description: >-
  Step by step tutorial with detail explanations
author: leo_choee
date: 2025-06-27
categories: [Github, Tutorial]
tags: [getting started]
pin: false
media_subpath: '/posts/20180809'
math: true
mermaid: true
---

> Reference book: Understanding Delta-Sigma Data Converters 2nd ed.

This page includes
1. Modifying high-pass filter to be suited for noise transfer function(NTF) in delta-sigma modulator (DSM)
2. Optimizing the selected NTF for better SQNR performance

## Step 1. Start with any high-pass filter you want

In this tutorial, we’ll use the Butterworth filter provided by MATLAB. For a deeper understanding, it’s highly recommended to refer to Chapter 4 of [Understanding Delta-Sigma Data Converters].

```matlab
[b, a] = butter(3,1/8,'high')
sys = tf(b,a,1,'Variable','z^-1')
```

$$
H\left(z\right)=\frac{0\ldotp 6735-0\ldotp 0204z^{-1} +2\ldotp 0204z^{-2} -0\ldotp 6735z^{-3} }{1-2\ldotp 2192z^{-1} +1\ldotp 7151z^{-2} -0\ldotp 4535z^{-3} },\omega =\left(0,\pi \right) 
$$

Since `the passband gain of the filter shoud be 1`, $H\left(z=\infty \right)=1$, $H\left(z\right)$ has to be scaled by 1/0.6735.

```matlab
b = b/0.6735;
sys = tf(b,a,1,'Variable','z^-1') 
```

$$
H\left(z\right)=\frac{1-3z^{-1} +3z^{-2} -z^{-3}=\left(1-z^{-1}\right)^3 }{1-2\ldotp 2192z^{-1} +1\ldotp 7151z^{-2} -0\ldotp 4535z^{-3} }
$$

Now, let's confirm whether our filter is stable using **pole-zero map**.

```matlab
zplane(b,a) 
```

You can find that the poles are inside the unit circle, and on the right half plane, which means that our filter is stable. While the filter is ready for immediate use, `two key factors` remain that directly impact its performance.

1. OBG (out-of-band gain) when $\omega =\pi \;\left(z=-1\right)$
2. IBG (in-band gain) when $\omega =0\;\left(z=1\right)$

When we calculate IBG, we approximate $1-z^{-1} =2\sin \left(\frac{\omega }{2}\right)\approx \omega$ with low-frequency input.

```matlab
OBG = polyval(b,-1)/polyval(a,-1);
fprintf("OBG = %.2f\n", OBG);
IBG = 1/polyval(a,1);
fprintf("IBG = %.2f\n", IBG);
```

One might ask, ***"Isn't an IBG of 23.58 too high? Couldn't it even increase the in-band noise?"***
However, this overlooks a key point: the very small value  appears in the numerator, significantly reducing its effect within the signal band. Moreover, **the total integrated noise power remains unchanged—it's simply reshaped by the high-pass nature of the NTF.** In other words, our NTF is still doing its job!
The resulting SQNR can be obatined as below equation and it's code implementation.

$$
P_{q\ldotp \mathrm{noise}} =\frac{\Delta^2 }{12}\cdot \frac{1}{2\pi }\int_{-\frac{\pi }{\mathrm{OSR}}}^{\frac{\pi }{\mathrm{OSR}}} {\left|\mathrm{NTF}\left(e^{j\omega } \right)\right|}^2 d\omega
$$

$$
P_{\mathrm{signal}} =\frac{A^2 }{2}
$$

$$
\mathrm{SQNR}=10\log_{10} \left(\frac{P_{\mathrm{signal}} }{P_{q\ldotp \mathrm{noise}} }\right)
$$

```matlab
vfs = 2; % assume that VDD = 1V
OSR = 64;
N = 2^14;
quantizer_bit = 4;

lsb = vfs/(2^quantizer_bit);
w = linspace(-pi/OSR, pi/OSR, N);

NTF = polyval(b,z.^-1)./polyval(a,z.^-1);
NTF_mag = abs(NTF);

P_signal = ((vfs/2)^2)/2;
P_qnoise = (lsb^2 / 12) * (1 / (2*pi)) * trapz(w, NTF_mag.^2);

sqnr = 10*log10(P_signal/P_qnoise);

fprintf('SQNR = %.2f dB\n', sqnr);
```

So, how does our current SQNR look? Is it sufficient for the target application?
While the current value may seem acceptable, there’s still room for improvement—so let’s optimize the performance a bit more in Step 2.

## Step 2. Optimize the NTF

The only way to improve the filter's performance with a fixed topology is to relocate the poles and zeros. Let's start from changing the poles.

### Relocating the poles

The filter's performance can be improved by ***shifting the 3-dB corner further away from the signal band.*** The cut-off frequency is modified when the poles are changed, and is shifted further away when the poles are moved closer to the high-frequency range $\left(\omega =\pi \right)$.
Thankfully, since we have a ready-made filter, we can achieve our goal simply by adjusting its parameters.

```matlab
% 3-dB corner from pi/8 to pi/4
[b, a] = butter(3,1/4,'high'); 
sys = tf(b,a,1,'Variable','z^-1')
```

Follow the procedure above in `Step 1`. You can use the prepared functions in [my git repository]

```matlab
b = b/0.4459;
sys = tf(b,a,1,'Variable','z^-1')
% zplane(b,a)
OBG = polyval(b,-1)/polyval(a,-1);
fprintf("OBG = %.2f\n", OBG);
IBG = 1/polyval(a,1);
fprintf("IBG = %.2f\n", IBG);
sqnr = calculate_NTF_SQNR(b, a, vfs, OSR, N, quantizer_bit); % function in git repository
fprintf('SQNR = %.2f dB\n', sqnr);
```

The SQNR is increased as 118.97 dB, which is about 15 dB higher than before.

This is a good point to **verify the rough location of our filter's cut-off frequency.** We can estimate the cut-off frequency by calculating the angle of the pole that is closest to the DC point (1, 0) on the z-plane.

```matlab
poles = pole(sys)
omega = atan(real(poles(1))/imag(poles(1)));
fprintf("ω = %.4f (derived)\nω = %.4f (𝜋/4)", omega, pi/4); % unit: rad/sample
```

Considering that the NTF consists of multiple poles and zeros, the discrepancy between the two values is small enough to conclude that the cut-off frequency is appropriately located as intended.

We can also intuitively compare both the OBG and IBG to their previous values by calculating the distances from the poles to (−1,0) and (1,0), respectively.

```matlab
OBG_rough = (1-(-1))^3/(-1-real(poles(1)))^2+(-imag(poles(1))^2);
IBG_rough = 1/(1-real(poles(1)))^2+(-imag(poles(1))^2);
OBG_rough_pre = (1-(-1))^3/(-1-0.7755)^2+(-0.2782^2);
IBG_rough_pre = 1/(1-0.7755)^2+(-0.2782^2);
fprintf("Previous: IBG = %.2f / OBG = %.2f \nOptimized: IBG = %.2f / OBG = %.2f" ...
    , IBG_rough_pre, OBG_rough_pre, IBG_rough, OBG_rough);
```

The result shows that the IBG is decreased, and the OBG is increased, as intended. These changes are also easily observable in the pole-zero plot.

### Spreading the zeros

***Will be updated soon!***

[Understanding Delta-Sigma Data Converters]: https://ieeexplore.ieee.org/servlet/opac?bknumber=5264508
[my git repository]: https://github.com/Leo-Choee