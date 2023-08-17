---
layout: post
title: a post with table of contents
date: 2023-03-20 11:59:00-0400
description: an example of a blog post with table of contents
categories: sample-posts toc
giscus_comments: true
related_posts: false
toc:
  beginning: true
---
This post shows how to add a table of contents in the beginning of the post.

# MP-Gadget Code Notes
Useful resources for understanding cosmological simulations:

- [Alexander Knebe's lecture notes](http://popia.ft.uam.es/aknebe/page3/files/ComputationalCosmology/)
- [gadget4 paper](https://wwwmpa.mpa-garching.mpg.de/gadget4/gadget4-code-paper.pdf)
- [Quinn97: Time stepping N-body simulations](https://arxiv.org/pdf/astro-ph/9710043.pdf)
- [Bagla04: Cosmological N-Body Simulations](https://arxiv.org/pdf/astro-ph/0411730.pdf)
- [Dolag08: Simulation techniques for cosmological simulations](https://arxiv.org/pdf/0801.1023.pdf)
- [Yepes Nbody Slides and references therein](https://workshops.ift.uam-csic.es/uploads/charla/275/Yepes_nbody.pdf)
- [Knebe's lecture](http://popia.ft.uam.es/ACO/links.html)

## Overview
(simplified pseudo-code for the main loop in `run.c`)

```
-------- Initializations ------------
hci_init (what is HCI?)
get_unit_system (independent of a)
gravpm_init_periodic
DriftKickTimes times = init_driftkicktime(ti_init)
atime = get_atime(times.Ti_Current)

----------- Main Loop -------------
times.Ti_Current = find_next_kick(times.Ti_Current, times.mintimebin);
atime = get_atime(times.Ti_Current);
if pm and all.randomoffset:
	update_random_offset(PartManager, rel_random_shift, All.RandomParticleOffset)

#--------- Drift/Ddcomp ----------
if PM_step or active_timebin:
	drift_all_particles(Ti_Last, times.Ti_Current, &All.CP, rel_random_shift)
	domain_decompose_full(ddecomp)
	
update_lastactive_drift(&times);

build_active_particles(&Act, &times, NumCurrentTiStep, atime);

#----------- Hydro ---------------
density(&Act, 1, DensityIndependentSphOn(), All.BlackHoleOn, times, &All.CP, &sph_predicted, GradRho_mag, &gasTree)
# update smoothing lengths in tree
force_update_hmax(Act.ActiveParticle, Act.NumActiveParticle, &gasTree, ddecomp);
hydro_force(&Act, atime, &sph_predicted, times, &All.CP, &gasTree);

#---------- PM --------------------
force_tree_full(&Tree, ddecomp, HybridNuTracer, All.OutputDir);
gravpm_force(&pm, &Tree, &All.CP, atime, units.UnitLength_in_cm, All.OutputDir, header->TimeIC, All.FastParticleType);

#--------- Tree ------------------
force_tree_full(&Tree, ddecomp, HybridNuTracer, All.OutputDir);
grav_short_tree(&Act, &pm, &Tree, NULL, rho0, HybridNuTracer, All.FastParticleType, times.Ti_Current);

#--------- First Half Kick ------------------
apply_half_kick(&Act, &All.CP, &times, atime); # both grav and hydro
update_kick_times(&times);
if pm:
	apply_PM_half_kick(&All.CP, &times);
	
#------- Cooling and extra physics --------
metal_return();
fof_seed(); # seed BH
do_heiii_reionization();
calculate_uvbg();
blackhole();
cooling_and_starformation();
force_tree_free(&Tree);

#--------- Output ----------------
if(WriteSnapshot):
	write_checkpoint();
	fof_save_groups();
 write_cpu_log();

#---------- Second half kick --------
find_timesteps(&Act, &times, atime, All.FastParticleType, &All.CP, asmth, NumCurrentTiStep == 0);
apply_half_kick(&Act, &All.CP, &times, atime);
update_kick_times(&times);
apply_PM_half_kick(&All.CP, &times);

```

## I/O
Driver: 

`checkpoint.c/write_checkpoint(int snapnum, int WriteGroupID, int MetalReturnOn, double Time, const Cosmology * CP, const char * OutputDir, const int OutputDebugFields)`

`petaio.c/petaio_save_snapshot(const char * fname, struct IOTable * IOTable, int verbose, const double atime, const Cosmology * CP)`





## Time integration

Detailed time integration scheme can be find in Section 4 of the [gadget4 paper](https://wwwmpa.mpa-garching.mpg.de/gadget4/gadget4-code-paper.pdf), also see [Quinn1997](https://arxiv.org/pdf/astro-ph/9710043.pdf):
$$
\mathbf{x}\left(\tau_n+\Delta \tau\right)=\mathbf{x}\left(\tau_n\right)+\mathbf{p}_n\left(\tau_n\right) \int_{\tau_n}^{\tau_n+\Delta \tau} \frac{\mathrm{d} \tau}{a^2 H(a)}
$$

$$
\mathbf{p}\left(\tau_n+\Delta \tau\right)=\mathbf{p}\left(\tau_n\right)+\mathbf{a}_n\left(\tau_n\right) \int_{\tau_n}^{\tau_n+\Delta \tau} \frac{\mathrm{d} \tau}{a H(a)}
$$

where $\mathbf{p} = a^2\dot{\mathbf{x}}$, $\mathbf{a} = a\dot{\mathbf{p}}$, $\tau = log(a)$, $dt = dlog(a)/H(a)$.

**Relation with physical quantities**

$\mathbf{r},\mathbf{v},\mathbf{f}$ are physical distance, velocity and acceleration; $\mathbf{x},\mathbf{p},\mathbf{a}$ are code distance, velocity and acceleration.
$$\mathbf{r} = a \mathbf{x}$$
$$\mathbf{v} = d\mathbf{r}/dt = \frac{1}{a}\mathbf{p} + \dot{a}\mathbf{x}$$
$$\mathbf{f} = d\mathbf{v}/dt = \dot{\mathbf{p}}/a + \ddot{a}\mathbf{x} = \frac{1}{a^2} \mathbf{a} + \ddot{a} \mathbf{x}$$

**[ChaNGa Notes](https://github.com/N-BodyShop/changa/wiki/ChaNGa-User-Guide#changa-with-black-holes)**
> Here be Dragons.
> 
> Going from comoving distance to physical distance is straightforward: r_phys = a r_cm. Physical velocity requires including the Hubble flow: v_phys = a (v_cm + H r_cm), where H is the hubble constant at the given redshift. Beware that v_circular is calculated differently, and therefore scales differently: v_c_phys = v_c_cm/sqrt(a). Finally, sometimes you want to know the ratio of the density to the critical density (e.g. to calculate the virial radius). This depends on cosmology. For Lambda + Matter critical Universes, we have rho/rho_c = rho_cm/(a^3 Omega_Lambda + Omega_matter), where the Omegas are present values and rho_c at the present is 1 as described above.
	
**In MP-gadget**

`get_exact_draft_fac` called in `drift.c`,`exchange.c`,`hydra.c`,`lightcone.c`

`get_exact_gravkick_fac` called in `blackhole.c`,`density.c`,`hydra.c`,`timestep.c`,`wind.c`


In `drift.c`

	ddrift = get_exact_drift_factor(CP, ti0, ti1);
	pp->Pos[j] += pp->Vel[j] * ddrift + random_shift[j];
	
        
In `timefac.c`:

`drift_integ(a)` = $1/H(a)/a^3$
	
`gravkick_integ(a)` = $1/H(a)/a^2$
	
`hydrokick_integ(a)` = $1/H(a)/a^{3(\gamma-1)+1}$
	
`get_exact_factor(a)`: integrate draft/kick factor from t0 to t1
$$\int_{a0}^{a1} fac(a)\,da$$

	a0 = exp(loga_from_ti(t0))
	
#### note on hubble
`CP->HubbleParam`: little h unit of 100km/s/Mpc

`CP->Hubble`: `CP->Hubble = HUBBLE * units.UnitTime_in_s`, 100km/s/Mpc in internal unit; set to 0.1 by default


`HUBBLE`: `#define  HUBBLE 3.2407789e-18 /* in h/sec */` = 100 km/s/Mpc in cgs

`HUBBLE * CP->HubbleParam`: hubble constant at z=0 in cgs

`CP->Hubble * CP->HubbleParam`: hubble constant at z=0 in internal time unit (UnitTime_in_s, 978 Myr)

`hubble_func(a) * CP->HubbleParam`: hubble constant at z in internal time unit

`UnitTime_in_s = 3.08568e+16`: chose such that `CP->Hubble=0.1`?





## Timestepping

### Timeline/Conversions

**Basic variables**

The simulated timespan is mapped onto the integer interval [0,TIMEBASE]


*Each integer time stores in the first 10 bits the snapshot number; Then the rest of the bits are the standard integer timeline, which should be a power-of-two hierarchy*

timebin是说bin的size

```
TIMEBINS=20
TIMEBASE = 2^21
MAXSNAPSHOTS = 2^11
dti_from_timebin(bin) = 2^(bin)
```

**Sync points and output**

2 << 1 = 4 and 2 << 2 = 8; 4 >> 1 = 2 (and 5 >> 1 = 2 since you round down);
The bitwise representation of a is shifted left b bits. This is the same as multiplying by (2 to the power of b).

```
	SyncPoints[0].a = TimeIC;
	SyncPoints[0].loga = log(TimeIC);
	
	SyncPoints[NSyncPoints].a = TimeMax;
    SyncPoints[NSyncPoints].loga = log(TimeMax);
        
    for(i = 0; i < NSyncPoints; i++) {
        SyncPoints[i].ti = (i * 1L) << (TIMEBINS);
    }
```
`SyncPoints[i].ti`: first 11 bits=i; last 21 bits=0.

`setup_sync_points(Cosmology * CP, double TimeIC, double TimeMax, double no_snapshot_until_time, int SnapshotWithFOF)`: initialize all sync points, their time, loga, writeFOF, etc.. `no_snapshot_until_time` ensures that we do not write output before the restart snapshot, even if those redshifts appears in the output list.


**Integer timeline and loga**

*Basically, the timeline is chopped into chunks based by syncpoints (snapshots), then the timeline between two consecutive snapshots is further divided into 2^21 ticks, everything is then in units of these ticks.*

`Dloga_interval_ti(inttime_t ti)`: returns the dloga corresponding to dti=1 at the current ti (valid between last snap and next snap), i.e. `(SyncPoints[lastsnap+1].loga - SyncPoints[lastsnap].loga)/2^21`; this number will differ from snap to snap.

`loga_from_ti(inttime_t ti)`: returns the current loga, calculated from `SyncPoints[lastsnap].loga + dti * Dloga_interval_ti(ti)`, where `dti = ti - SyncPoints[lastsnap].ti` = the last 21 bits of ti?

`dloga_from_dti(inttime_t dti, const inttime_t Ti_Current)`: get the dloga from dti according to the current snapshot interval: `returns Dloga_interval_ti(Ti_Current) * dti` *but what if dti crosses the next snap?*

`get_dloga_for_bin(int timebin, const inttime_t Ti_Current)`:
returns the length of the time bin: `Dloga_interval_ti(Ti_Current) * 2^timebin`

`init_timebins(TimeInit) = ti_from_loga(log(TimeInit))`

`int is_timebin_active(int i, inttime_t current)` 如果current可以被第i个bin的长度整除（`dti_from_timebin(i)`），**timebin**而不是ti就active (ti=12的话，bin1，2 active，3不active). 
**用法**：`is_timebin_active(pp->TimeBinHydro, ti_current)`: 判断在某一个time ti，某一个particle的特定calculation是不是active.




### Timestep calculations
**KD times**

`find_timesteps` (called before the second half kick): assign new timesteps to the active particles, now that we know they have synched TiKick and TiDrift, and advance the PM timestep.

`apply_half_kick(const ActiveParticles * act, Cosmology * CP, DriftKickTimes * times, const double atime)`
对每一个bin （1-20），`new_kick = times->Ti_kick[bin] + dti_from_timebin(bin)/2`

For each bin, plan to kick from `times->Ti_kick[bin]` to `times->Ti_kick[bin] + dti_from_timebin(bin)/2`, the kick factors are stored in `gravkick[bin]`

```
if(is_timebin_active(bin_gravity, times->Ti_Current)) {
   		do_grav_short_range_kick(&P[i], P[i].FullTreeGravAccel, gravkick[bin_gravity]);
   		do_hydro_kick(i, dt_entr, gravkick[bin_hydro], hydrokick[bin_hydro], atime);
```
   
 Note: `apply_hydro_half_kick` done separately only under hierarchical gravity

**grav-timestep**

Basically a small fraction of softening/acc

`get_timestep_gravity_dloga`:  $dt = \sqrt{2 \epsilon_{\rm tol} \, a \, \epsilon_g / \mathbf{a}_{\rm grav}}$

**hydro-timestep**

*ref: [gadget4 paper](https://wwwmpa.mpa-garching.mpg.de/gadget4/gadget4-code-paper.pdf) Eq. 77-81*

- local particle-based Courant-Friedrichs-Lewy (CFL) hydrodynamic timestep (smoothing length/sound speed)

$$
\Delta t_i^{\mathrm{cfl}}=C_{\mathrm{CFL}} 2 h_i / v_i^{\mathrm{sig}, \max }
$$

- timestep criterion that restricts the al- lowed rate of change of the smoothing length (or equivalently den- sity) per timestep

$$
\Delta t_i^{\mathrm{dens}}=C_{\mathrm{CFL}} \frac{h_i}{\left|\mathrm{~d} h_i / \mathrm{d} t\right|},
$$


**long-range timestep**
`get_long_range_timestep_dloga`: just set PM step to maxtimebin for non-cosmological sim; pm step will be shortened to the longest tree timestep in find_timestep().






## Gravity


Basic equation of motion (x: comoving) ([reference1](https://arxiv.org/pdf/astro-ph/0411730.pdf),[reference2](https://arxiv.org/pdf/0801.1023.pdf)):


$$\ddot{x} + 2H(a)\dot{x} = \frac{1}{a}\,\nabla\phi$$

$$\nabla^2\phi = 4\pi G \bar{\rho}(t)a^2\delta = \frac{3}{2}H_0 ^2\Omega_0 \frac{\delta}{a}$$

$$
\frac{\mathrm{d} \mathbf{v}}{\mathrm{d} t} =-\frac{\nabla \Phi}{a}-\mathbf{v} H(a)
$$

where $\mathbf{v}$ is proper peculiar velocity.

### PP Force


In `gravshort_tree.c`:

`apply_accn_to_output`:

	fac = mass / (r2 * r);
	for(i = 0; i < 3; i++)
        output->Acc[i] += dx[i] * fac;
        output->Potential += facpot;
        

$$
\phi(\mathbf{x})= -\sum_{j=1}^N \frac{m_j}{\left|\mathbf{x}_j-\mathbf{x}+\mathbf{q}_j^{\star}\right|+\epsilon\left(\left|\mathbf{x}_j-\mathbf{x}+\mathbf{q}_j^{\star}\right|\right)} 
+\sum_{j=1}^N m_j \psi\left(\mathbf{x}_j-\mathbf{x}+\mathbf{q}_j^{\star}\right)
$$

	
### PM Force	

Reference: [SUPERBOX paper](https://core.ac.uk/download/pdf/25290081.pdf)

Possion's equation in real space (a convolution of density rho with the Green's function H):

$$
\Phi_{i j k}=G \sum_{a, b, c=0}^{N-1} \varrho_{a b c} \cdot H_{a-i, b-j, c-k}, \quad i, j, k=0, \ldots, N-1
$$

where the Green's function for gravivational potential is:

$$
\begin{aligned}
H_{i j k} &=\frac{1}{\sqrt{i^2+j^2+k^2}}, \quad i, j, k=0,1, . ., N \\
H_{000} &=4 / 3 .
\end{aligned}
$$




### units and scale factors
Scale factors

`P[i].Vel[k]` <- $a\,v_{\rm pec}$;

`P[i].Pos[k]` <- $x_{\rm prop}/a$;

`P[i].DFAccel[k]` <- $a^2\,a_{\rm df}$;

Code units

velocity: $a\,km/s = 10^3\, m/s\;\times a$

length: $kpc/h/a = 3.086\times 10^{19}\, m\;/h/a$;

time: length/velocity = $[kpc]/[km/s] = 3.086e16\,s = 978\,Myr$; (any extra a or h?)

acc: $a^2\times length/time^2 = [km/s]^2/[kpc] = 3.24 \times 10^{-14}\,m/s^2\;\times a^2$

rho: $mass/length^2 = 1e10/h\,M_\odot/[kpc/h/a]^3 = 10^{10}\,a^3\,h^2 M_\odot/kpc^3 = 10\,a^3\,h^2 M_\odot/pc^3 = a^3\,h^2\;6.77\times 10^{-19}\,kg/m^3$



## Hydro/SPH
**basic thermo relations**

Equation of state
$$
P_i=u_i(\gamma-1) \rho_i
$$

Entropy function A
$$
P_i=A_i(s) \rho_i^\gamma
$$

Energy, entropy, density, and pressure
$$
u \equiv \bar{P}_i /\left[(\gamma-1) \bar{\rho}_i\right]=(\gamma-1)^{-1} A_i \bar{\rho}_i^{\gamma-1}
$$

**Euler equations** (see [here](https://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/Chapter_12.pdf))

$$
\begin{gathered}
\frac{\mathrm{d} \rho}{\mathrm{d} t}+\rho \nabla \cdot \mathbf{v}=0, \\
\frac{\mathrm{d} \mathbf{v}}{\mathrm{d} t}+\frac{\nabla P}{\rho}=0, \\
\frac{\mathrm{d} u}{\mathrm{~d} t}+\frac{P}{\rho} \nabla \cdot \mathbf{v}=0,
\end{gathered}
$$

where $\mathrm{d} / \mathrm{d} t=\partial / \partial t+\mathbf{v} \cdot \nabla$ is the convective derivative.

**SPH formulation**

Pressure-entropy formulation
$$
\bar{P}_i=y_i^\gamma=\left[\sum_{j=1}^N m_j A_j^{1 / \gamma} W_{i j}\left(h_i\right)\right]^\gamma
$$

sound speed in pressure-entropy sph:
$$
c_i = \sqrt{\gamma\,P_i/\rho_{u,i}}
$$
standard sph
$$
c_i = \sqrt{\gamma\,P_i/\rho_i}
$$

EOM:
$$
\begin{aligned}
\frac{\mathrm{d} \mathbf{v}_i}{\mathrm{~d} t} & =-\sum_{j=1}^N m_j\left(A_i A_j\right)^{\frac{1}{\gamma}}\left[\frac{f_i \bar{P}_i}{\bar{P}_i^{2 / \gamma}} \nabla_i W_{i j}\left(h_i\right)+\frac{f_j \bar{P}_j}{\bar{P}_j^{2 / \gamma}} \nabla_i W_{i j}\left(h_j\right)\right] \\
f_i & =\left[1+\frac{h_i}{3 \bar{P}_i^{1 / \gamma}} \frac{\partial \bar{P}_i^{1 / \gamma}}{\partial h_i}\right]^{-1}
\end{aligned}
$$

### entropy
in `init.c`: `SphP[i].Entropy = GAMMA_MINUS1 * u_init / pow(SphP[i].EgyWtDensity / a3 , GAMMA_MINUS1);`

in 	`winds.c`: `SPHP(other).Entropy += therm/enttou;`

in `sfr_eff.c`: `enttou = entropy_to_u(SPHP(i).Density, a3inv);`

in `blackhole.c`: `enttou = pow(SPHP(other).Density * BH_GET_PRIV(lv->tw)->a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;`

$$s2u = (\rho/a^3)^{\gamma-1}/(\gamma - 1)$$ (same as the basic thermal relation above)

`O->EgyRho += mass_j * EntVarPred * wk; SPHP(i).EgyWtDensity /= EntPred;`


### density estimation (density.c)
Initialize the kernel by type; calculate the gas particles within the kernel (depends on kernel type, 113 for quintic)

```double
GetNumNgb(enum DensityKernelType KernelType)
{
    DensityKernel kernel;
    density_kernel_init(&kernel, 1.0, KernelType);
    return density_kernel_desnumngb(&kernel, DensityParams.DensityResolutionEta);
}
```

`EntVarPred = exp(1./GAMMA * log(EntVarPred));` why??




### hydro(hydra.c)
`PressurePred(MyFloat EOMDensityPred, double EntVarPred)`: `return pow(EntVarPred * EOMDensityPred, GAMMA);` 
This is likely $\bar{P}_i=y_i^\gamma=\left[\sum_{j=1}^N m_j A_j^{1 / \gamma} W_{i j}\left(h_i\right)\right]^\gamma$

Then ` HYDRA_GET_PRIV(tw)->PressurePred[i] = PressurePred(SPH_EOMDensity(&SphP[i]), SPH_predicted->EntVarPred[i]);`
`SPH_EOMDensity` returns either Density or EgyWtDensity

Key:
`hfc += P[other].Mass *
            (dwk_i*iter->p_over_rho2_i*EntVarPred/I->EntVarPred +
            dwk_j*p_over_rho2_j*I->EntVarPred/EntVarPred) / r;`

hydro force added in `hydro.c`: 
`O->Acc[d] += (-hfc * dist[d]);` dist: distance

`TREEWALK_REDUCE(SPHP(place).HydroAccel[k], result->Acc[k]);`

### 




## Heating, Cooling, ionization equilibrium
In cooling.c

```
double
get_neutral_fraction_phys_cgs(double density, double ienergy, double helium, const struct UVBG * uvbg, double * ne_init)                                                              
{
    double logt;
    double ne = get_equilib_ne(density, ienergy, helium, &logt, uvbg, *ne_init);
    double nh = density * (1-helium);
    double photofac = self_shield_corr(nh, logt, uvbg->self_shield_dens);
    *ne_init = ne/nh;
    return nH0_internal(logt, ne, uvbg, photofac);
```

Set of ionization equilibrium equations from KWH95:

$$
\begin{aligned}
n_{\mathrm{H}_0} & =n_{\mathrm{H}} \alpha_{\mathrm{H}_{+}} /\left(\alpha_{\mathrm{H}_{+}}+\Gamma_{\mathrm{eH}_0}+\Gamma_{\gamma \mathrm{H}_0} / n_{\mathrm{e}}\right) \\
n_{\mathrm{H}_{+}} & =n_{\mathrm{H}}-n_{\mathrm{H}_0}, \\
n_{\mathrm{He}_{+}} & =y n_{\mathrm{H}} /\left[1+\left(\alpha_{\mathrm{He}_{+}}+\alpha_{\mathrm{d}}\right) /\left(\Gamma_{\mathrm{eHe}_0}+\Gamma_{\gamma \mathrm{He}_0} / n_{\mathrm{e}}\right)+\left(\Gamma_{\mathrm{eHe}_{+}}+\Gamma_{\gamma \mathrm{He}_{+}} / n_{\mathrm{e}}\right) / \alpha_{\mathrm{He}_{++}}\right] \\
n_{\mathrm{He}_0} & =n_{\mathrm{He}_{+}}\left(\alpha_{\mathrm{He}_{+}}+\alpha_{\mathrm{d}}\right) /\left(\Gamma_{\mathrm{eHe}_0}+\Gamma_{\gamma \mathrm{He}_0} / n_{\mathrm{e}}\right) \\
n_{\mathrm{He}_{++}} & =n_{\mathrm{He}_{+}}\left(\Gamma_{\mathrm{eHe}_{+}}+\Gamma_{\gamma \mathrm{He}_{+}} / n_{\mathrm{e}}\right) / \alpha_{\mathrm{He}_{++}} \\
n_{\mathrm{e}} & =n_{\mathrm{H}_{+}}+n_{\mathrm{He}_{+}}+2 n_{\mathrm{He}_{++}} .
\end{aligned}
$$




## Star-formation
> "The star formation model is unchanged from Feng et al.
(2016a), itself an implementation of Springel & Hernquist
(2003). Gas is allowed to cool both radiatively (Katz et al.
1996) and via metal line cooling. We include self-shielding of
dense gas via the fitting function of Rahmati et al. (2013). We
approximate the metal cooling rate by scaling a solar metallicity template according to the metallicity of gas particles,
following Vogelsberger et al. (2014). We include a correction
for the formation of molecular hydrogen, and its effect on star
formation at low metallicities, according to the prescription
by Krumholz & Gnedin (2011). Stars are formed with 1/4 of
the mass of a gas particle. Since we have also implemented
mass return, the mass of a gas particle after forming 4 stars
may exceed zero. To prevent runaway enrichment as gas cycles between a star and gas particle, the 5th star spawned
from a gas particle contains the entire mass of the particle."




## winds, stellar feedbacks

> "A stellar wind feedback model (Okamoto et al. 2010) is included, which assumes wind speeds proportional to the local
one dimensional dark matter velocity dispersion σDM:
vw = κwσDM , (3)
where vw is the wind speed. κw is a dimensionless parameter, which we take to be 3.7 following Vogelsberger et al.
(2013). Winds are sourced by newly formed star particles,
which randomly pick gas particles from within their SPH
smoothing length to become wind particles. The total mass
loading is (vw/350km/s)−2 where 350 km/s is in physical
units. Once a particle is in the wind, it is hydrodynamically"

## Metal return

> "We include a model for the return of mass and metal to
the inter-stellar medium from massive stars. The general
approach follows Vogelsberger et al. (2013); Pillepich et al.
(2018a), although the implementation is independent and we
have created our own metal yield tables. Each star particle is
treated as a single stellar population, with the proportion of
massive stars set by a Chabrier initial mass function (IMF)
(Chabrier 2003). Massive stars have a finite lifetime, set using
the lifetime tables of Portinari et al. (1998), and at the end of
this lifetime they return mass and metal to neighbouring gas
particles proportional to their fraction of the underlying stellar population. Nine metal species are followed: H, He, C, N,
O, Ne, Mg, Si, Fe. The total metallicity is tracked separately.
Metal and mass yields from AGB stars, SnII and Sn1A are
pre-compiled and included as a function of mass and stellar
metallicity. Stars with masses 1M − 8M are assumed to
return via an Asymptotic Giant Branch (AGB) stage, with"





## Parallelization
Useful resources

- [ChaNGa MPI paper](https://charm.cs.illinois.edu/newPapers/08-03/paper.pdf)
