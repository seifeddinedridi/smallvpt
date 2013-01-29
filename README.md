smallvpt
========

<b>smallvpt</b> is a small volumetric path tracer written in 150 lines of code, it is an extension of <a href="http://www.kevinbeason.com/smallpt/">smallpt</a>, which is a tiny path tracer written by Kevin Beason.

It currently includes basic volumetric light transport supporting multiple scattering and homogeneous media. The medium is defined as a sphere surrounding the scene, with parameters <b>sigma_s</b> and <b>sigma_a</b> to control the scattering and absorption of light.

Feel free to play with it as you like, if you render something cool or spot an error, please let me know ;)

I apologize in advance for the noise; my computer is not fast enough to clear it, and I can't afford an overnight render due to the untolerable sound of a choking CPU :)

Without further ado, here are a couple of renders showing what it is capable of:

![My image](https://raw.github.com/D-POWER/smallvpt/master/Renders/Foggy%20Cornell%20Box%20-%20%5B10000spp%5D.png)

![My image](https://raw.github.com/D-POWER/smallvpt/master/Renders/Volumetric%20caustics%20-%20%5B10000spp%5D.png)

![My image](https://raw.github.com/D-POWER/smallvpt/master/Renders/Glass%20in%20a%20medium%20-%20%5B1024spp%5D.png)

![My image](https://raw.github.com/D-POWER/smallvpt/master/Renders/image%20-%200.01%20sigma_s%20%5B2048spp%5D.png)
