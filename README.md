# Foregrounds

Code for generating a box with foreground contamination.

## The model
This model is for extra-galactic foregrounds. The number of foregrounds per pixel is specified by the `foregrounds_per_pix` flag, where the default is 5. The temperature spectrum of every foreground is given by T(nu) = A*(nu/150)**(-alpha), and every foreground has a different spectrum. For each foreground: 

- A is drawn from a power law distribution of the form P(x) = N*x**(-a), where the index a is specified, and the defaults is `a = 4.`. The mean of pixels at 150MHz is also specified by `mean_150MHz` and the default is 300 K. The minimum temperature of the distribution P(x) is determined to satisfy the index and mean.

- alpha is drawn from a gaussian distribution of specified mean and std. The defaults are `mean, std = 2.5, 0.5`

## The output
There are two different types of output boxes. 

### 1. Frequency space output

This gives a box that is uniformly spaced in frequency. We can create it the following way:
```python
f = Foregrounds()
foregrounds, frequencies = f.frequency_space(bandwidth = [nu_min,nu_max])
```

foregrounds contains a n by n by n array. The range of frequencies is given by the bandwidth and the spacing by the size of the box. The `frequencies` variable contains the considered frequencies. `foregrounds[i]` will contain a n by n an array of what the sky would look like when observed at a frequency `frequencies[i]`. For example, at 100 MHz, for one foreground per pixel: 

![alt text](https://i.imgur.com/ZeKYcBI.png)


### 3. Real space output

This gives a box that is uniformly spaced in real space. We can create it the following way:
```python
f = Foregrounds()
foregrounds = f.real_space(central_freq = nu_0, L = L)
```
`foregrounds` then contains a box centered on the comoving distance that corresponds to `nu_0`, of side length specified by `L`. 
