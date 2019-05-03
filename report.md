# Computer Graphics Report

## Raytracer

##### Basic version
<p float="left">
  <img src="raytracer/screenshots/basic.png" width="450" />
</p>

##### Optimisations
- Cramer's rule for computing triangle intersections
- OpenMP parallelism

##### Phong light
- Glossy surfaces are modeled using Phong light
<p float="left">
  <img src="raytracer/screenshots/phong_specular.png" width="260" />
  <img src="raytracer/screenshots/phong_glossy.png" width="260" />
</p>

##### Custom models
<p float="left">
  <img src="raytracer/screenshots/teapot.png" width="260" />
  <img src="raytracer/screenshots/shuttle.png" width="260" />
</p>
<p float="left">
  <img src="raytracer/screenshots/pine_tree.png" width="172" />
  <img src="raytracer/screenshots/trumpet.png" width="172" />
  <img src="raytracer/screenshots/ursache.png" width="172" />
</p>

#### Monte Carlo Path Tracing
- Point light
- Anti-aliasing and soft shadows by shooting slightly perturbed rays each time
<p float="left">
  <img src="raytracer/screenshots/point_light_400spp.png" width="450" />
</p>

##### Emissive materials
- Point lights are replaced by emissive surfaces   
<p float="left">
  <img src="raytracer/screenshots/hdr_14kspp.png" width="450" />
</p>

##### Tone mapping
- The color values can sometimes be higher than `vec3(1,1,1)`, which means they are clipped when passed to `PutPixelSDL`. To convert the colors from HDR to LDR we experimented with different tone mapping and gamma correction methods.
<p float="left">
  <img src="raytracer/screenshots/final3_25kspp.png" width="450" />
</p>

##### BRDF for specular (mirrors) and glossy materials
<p float="left">
  <img src="raytracer/screenshots/mirror_5kspp.png" width="260" />
  <img src="raytracer/screenshots/glossy_2kspp.png" width="260" />
</p>

##### Transparent & refractive materials
- The amount of light reflected when light enters another medium is computed using Schlick's approximation for Fresnel's equations
- Total internal refraction
<p float="left">
  <img src="raytracer/screenshots/refraction_5kspp.png" width="450" />
</p>

##### Beer-Lambert law
- Beer’s law models light being absorbed over distance as it travels through a transparent material. The extinction (absorption) describes how much of each colour channel absorbs over distance. For example, a value of `(8.0, 2.0, 0.1)` would make red get absorbed the fastest, then green, then blue, so would result in a blueish, slightly green color object. The first picture shows a fully transparent object, while the second cube has a refraction index of `1.52`.
<p float="left">
  <img src="raytracer/screenshots/beer_5kspp.png" width="450" />
</p>
<p float="left">
  <img src="raytracer/screenshots/blue_glass_10kspp.png" width="450" />
</p>

##### Fog density
-  Based on the distance travelled by a ray and on the height
<p float="left">
  <img src="raytracer/screenshots/fog_low_18kspp.png" width="450" />
</p>
<p float="left">
  <img src="raytracer/screenshots/fog_high_7kspp.png" width="450" />
</p>

#### Final renders
![basic](raytracer/screenshots/final3_25kspp.png)
![basic](raytracer/screenshots/final1_50kspp.png)
![basic](raytracer/screenshots/final2_20kspp.png)

## Rasteriser

##### Object loader
<p float="left">
  <img src="rasteriser/screenshots/lucy_28mil.png" width="450" />
</p>
<p float="left">
  <img src="rasteriser/screenshots/ursache.png" width="172" />
  <img src="rasteriser/screenshots/al.png" width="172" />
  <img src="rasteriser/screenshots/cessna.png" width="172" />
</p>
<p float="left">
  <img src="rasteriser/screenshots/pacanele.png" width="172" />
  <img src="rasteriser/screenshots/lamp.png" width="172" />
  <img src="rasteriser/screenshots/skyscraper.png" width="172" />
</p>

##### Texture mapping
- Used barycentric coordinates
![basic](rasteriser/screenshots/screenshot_spot_alb.png)

##### FXAA anti-aliasing
- Box without AA
![basic](rasteriser/screenshots/wo_fxaa.png)
- Box with AA
![basic](rasteriser/screenshots/w_fxaa.png)
