# Foregrounds

Code for generating a box with foreground contamination.

## The model
This model is for extra-galactic foregrounds. The number of foregrounds per pixel is specified by the `foregrounds_per_pix` flag, where the default is 5. The temperature spectrum of every foreground is given by T(nu) = A*nu**(-alpha), and every foreground has a different spectrum. For each foreground: 

- A is drawn from a power law distribution of the form P(x) = N*x**(-a), where the minimum temperature and the index a are specified. The defaults are `tmin, a = 3., 4.`. 

- alpha is drawn from a gaussian distribution of specified mean and std. The defaults are `mean, std = 2.5, 0.5`

## The output
To access the temperatures: 

```python
f = Foregrounds()
foregrounds = f.frequency_space_temps
```

foregrounds contains a n by n by n array. The range of frequencies is given by f.nu, and foregrounds[i] will contain a n by n an array of what the sky would look like when observed at a frequency f.nu[i]. For example, at 100 MHz: 

![alt text](https://i.imgur.com/ZeKYcBI.png)
