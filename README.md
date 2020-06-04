# ACSE-4-armageddon

Asteroids entering Earth’s atmosphere are subject to extreme drag forces that decelerate, heat and disrupt the space rocks. The fate of an asteroid is a complex function of its initial mass, speed, trajectory angle and internal strength. 

[Asteroids](https://en.wikipedia.org/wiki/Asteroid) 10-100 m in diameter can penetrate deep into Earth’s atmosphere and disrupt catastrophically, generating an atmospheric disturbance ([airburst](https://en.wikipedia.org/wiki/Air_burst)) that can cause [damage on the ground](https://www.youtube.com/watch?v=tq02C_3FvFo). Such an event occurred over the city of [Chelyabinsk](https://en.wikipedia.org/wiki/Chelyabinsk_meteor) in Russia, in 2013, releasing energy equivalent to about 520 [kilotons of TNT](https://en.wikipedia.org/wiki/TNT_equivalent) (1 kt TNT is equivalent to 4.184e12 J), and injuring thousands of people ([Popova et al., 2013](http://doi.org/10.1126/science.1242642); [Brown et al., 2013](http://doi.org/10.1038/nature12741)). An even larger event occurred over [Tunguska](https://en.wikipedia.org/wiki/Tunguska_event), an unpopulated area in Siberia, in 1908. 

This tool predicts the fate of asteroids entering Earth’s atmosphere for the purposes of hazard assessment.
This tool is based on a numerical solver that uses the improved euler technique to solve a system of differencial equations in order to calculate the fate.

While this tool is fast, once an asteroid has passed the atmosphere there would be no time to evacuate the area that will be affected in the timescales involved. This is a job for the astronomers, since it would take days to evacuate any large areas. 

The real value of this code, is to help us determine __what exactly happened__ and cooridnate the emergency services. It would be expected that for a catastrophic impact, communications would be severed. With this tool we would know exactly the area impacted (x) and we will have a good estimate of the damage, since we know __the energy released__ and __where__ the asteroid burst.

### Functionality


**Chelyabinsk**
We have used our solver to perform an inversion on the data from the 2014 Chelyabinsk asteroid incident. 
We had a list of altitudes and energy loss so the first oreder of bussiness was to convert those descrete data points to a continious function in order to compare the results more accurately. We performed this using the scipy-intergrate function. An example is seen from the figure bellow:
https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/images/approximation%20of%20cher.png

We inverted by running our solver with multiple combinations of Strenght and Radious as input data. We took the energy that resulted and calculates its RMS compared with the raw data. We plotted the results bellow:
https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/images/inversion.png

https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/images/inversion_zoomed.png

Unfortunately the numerical solver was not very accurate, so the curve was way off from the chelyabinsk curve. The least error was given with a Radious of 10m and strength of 10000000.


**Ensemble of Parameter spaces**
Arguments given by the user are being varied in the probability analysis. Then, for the given parameters program calculates what is the probability that a asteroid with given parameters will appear and then calculates the altitude on which the asteroide will burst. 
The results on the graph below are for the analysis with varied density of the asteroid. The density is being varied from 1000 kg/m^3 to 7000 kg/m^3 with 1000 data points. For each data point a probability of such density occuring is being calculated. Then, the program calculate the burst altitude for a given combination of parameters (though in this case everything but density remain unchanged).
https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/images/4000_density.png

The second plot presents an analysis where radius is varied from 5 to 15m, angle is being varied from 15 to 75 degrees and strength is varied from 10^3 to 10^7 Pa.

https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/images/15_radius_angle_strength.png


**Error Analysis**
In order to evaluate the numerical solver we compared it to the `scipy odeint` results. Our code gave near identical results with the scily solution as seen from the graph bellow:
https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/images/Image%20from%20iOS.png

**Accuracy - Tolerance Analysis**
This part compares 2 scenarios' (cratering & airburst) the analytical solution and the numerical solution with the same constant & parameters. The final plot shows the relationship between the tolerence and accuracy, which can tell an appropriate tolerence under chosen threshold accuracy. This is summarised in the figure bellow:
https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/images/image.png
### User instructions & Installation Guide

We have stived to make this tool as accessible as possible, therefore we have created a Graphical User Interface to go along with the code.

Simply install [git](https://www.linode.com/docs/development/version-control/how-to-install-git-on-linux-mac-and-windows/) on your computer.
Then open the command line and type:
    git clone https://github.com/acse-2019/acse-4-armageddon-hygiea.git
This will copy the reposity(folders) and then you need to open the *gui ipython* notebook. There you will be promted to enter the initial values of the asteroid you observed and the program will return back all the values for *Velocity, mass, energy, radious, lateral location, angle and altitude*. It will also print graphs for these parameters.


### Documentation

We have created a pdf of the documentation from sphinx and included that in the folder:
https://github.com/acse-2019/acse-4-armageddon-hygiea/blob/user_guide/acse42-documentation.pdf


### Testing

The tool includes several tests, which you can use to check its operation on your system. With [pytest](https://doc.pytest.org/en/latest) installed, these can be run with

```
python -m pytest armageddon
```
