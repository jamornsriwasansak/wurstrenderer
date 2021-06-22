# wurstrenderer (as of 2019)

![t2](https://github.com/jamornsriwasansak/wurstrenderer/blob/master/readme/tonemapped.jpg)

![t1](https://github.com/jamornsriwasansak/wurstrenderer/blob/master/readme/tea.jpg)

WurstRenderer is a CPU-based rendering system I made in early 2019. It focuses on experimentability and readability. For instance, the length of the [bidirectional path tracer code](https://github.com/jamornsriwasansak/wurstrenderer/blob/master/src/renderer/bidirpath.h) is only 430 lines long. All integrators (3d) are compared against reference generated from Mitsuba. The result of convergence plots indicates that the implementation of all integrators is correct.

> To avoid accidentally releasing code written by collaborators or sensitive code, all changes after mid 2019 were removed.

## Implemented:

In this framework, it supports:
* Real-time viewer
* Obj importing and partially PBRT importing.
* Oren-Nayar BSDF
* Mixed BRDF
* Modified Phong [[Lafortune and Willems 1994](http://www.lafortune.eu/publications/Phong.html)]
* Blue Dithering Mask Generation via Simulated Annealing [[Georgiev and Fajardo 2016](https://dl.acm.org/doi/abs/10.1145/2897839.2927430)]
* Bidirectional Path Tracing [[Veach's thesis 1997](https://graphics.stanford.edu/papers/veach_thesis)]
* Matrix Bidirectional Path Tracing [[Chaitanya et al. 2018](http://www.cim.mcgill.ca/~derek/files/mbdpt-e.pdf)] (haven't verified the correctness yet)
* Path Tracing Next Event Estimation with MIS [[Veach's thesis 1997](https://graphics.stanford.edu/papers/veach_thesis)]
* Direct RayTracing of Phong Tessellation [[Ogaki and Tokuyoshi 2011](http://www.jp.square-enix.com/tech/library/pdf/EGSR2011.pdf)]
* Precomputed Radiance Transfer (Diffuse) [[Sloan et al. 2002](https://sites.fas.harvard.edu/~cs278/papers/prt.pdf)]
* Precomputed Radiance Transfer (Glossy) [[Sloan et al. 2002](https://sites.fas.harvard.edu/~cs278/papers/prt.pdf)]
* Primary Sample Space Metropolis Light Transport [[Kelemen and Szirmay-Kalos
 2001](https://cg.informatik.uni-freiburg.de/intern/seminar/raytracing%20-%20Kelemen%20-%202002%20-%20Metropolis.pdf)]
* Environment Map with several representation includes equirectangular, spherical harmonics, spherical gaussians, and spherical fibonacci (thanks to inverse mapping [[Keinert et al. 2015](https://dl.acm.org/doi/10.1145/2816795.2818131)]).
* Spherical Gaussian Lobes Fitting via LibTorch
* Homogeneous Volumetric Rendering with Isotropic & Henyey-Greenstein phase functions.
> Note: It seems that homogeneous volumetric rendering is broken again...

## Flatland?:
Flatland is the world where two-dimensional creatures live in. 
It is quite useful for understanding complicated algorithms.
The renderer supports two integrators:
* Flatland Path Tracing
* Flatland PSSMLT

It simply takes the 3d scene and automatically take a slice along Z-Axis to create a 2d scene.

For flatland integrator, it can preview path density as well. I found it is very helpful for understanding MCMC & path guiding experiments.

![pt](https://github.com/jamornsriwasansak/wurstrenderer/blob/master/readme/pt_result.jpg) | ![bdpt](https://github.com/jamornsriwasansak/wurstrenderer/blob/master/readme/bdpt_result.jpg) | ![pt2d](https://github.com/jamornsriwasansak/wurstrenderer/blob/master/readme/pt2d_result.jpg)  | ![visualize](https://github.com/jamornsriwasansak/wurstrenderer/blob/master/readme/visualized.jpg)
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
Path Tracing NEE + MIS  | Bidirectional Path Tracing| Path Tracing in Flatland | Path Density Visualization

Scene format example:
```
{
    "viewer" : true,
    "camera" : 
    {
        "pos" : [4, 1, -2],
        "lookat": [-1, 1, 0],
        "up" : [0, 1, 0],
        "resolution" : [1024, 512],
        "fovy": 50
    },
    "meshes" : 
    [
        {
            "path" : "fireplace_room/fireplace_room.obj"
        }
    ],
    
    "render" : [
        {
            "integrator" : "nextevent_path_tracer",
            "output": "pt_result.pfm",
            "num_spp" : 100
        },
        {
            "integrator" : "bidirectional_path_tracer",
            "output": "bdpt_result.pfm"
        }
    ]
}
```

## Misc
It uses its own math library and automatically switch to SIMD if your CPU supports.
It strictly use double only. The performance is not its forte.
