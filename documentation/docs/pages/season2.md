# Season 2

## Navigating MOM6 as a non-oceanographer: a MOM6 mini-app for GPU work.

Date: 7/05/2026.

Presenter: Jorge Gálvez Vallejo (@JorgeG94).

Link to repo: [hackathon_mom6_miniapp](https://github.com/JorgeG94/hackathon_mom6_miniapp)

### Brief Intro

Joining a new project in a field you did not study can be extremely daunting! But fear not, for 
fortune favors the brave. I'd say there are three options here:

- You are an oceanographer that has little knowledge of high performance computing 
- You are a high performance computing person that has no knowledge of oceanography 
- You are a new student and don't know either

The third one naturally has the most to overcome, but they also have the most time. The first one
is simply "we have to learn how to code", the second one is "I have to learn oceanography". However, 
the second one might not need to learn a lot of oceanography since their contributions don't depend 
on them understanding the physics of the ocean dynamics. 

I am a stubborn person that doesn't like to write code I don't understand, so I basically like to 
create more work for myself.  

### Detailed notes (flow of consciousness) 

As an HPC person whose goal is to either optimize code or port it to different architectures the key 
aspect to understand is *what are the bottlenecks?* in the current code. In ocean simulations, to my 
annoyance, there is no _one_ routine that dominates 80% of the time. This makes my work harder, usually
one can quickly say "ah yes, the calculation of this physical thing takes most of the time" so one spends
the most time there. Then naturally you go ask the people that work on the codebase about the bottlenecks and
because everyone is using it for a different thing you will hear different things, this is where the oceanography
knowledge becomes useful! 

What is the best way to understand a gigantic codebase with thousands of lines of code? Simplify it! I.e. the creation of a mini-app. 
This is also the best thing to do for a hackathon, because you can capture algorithmic complexity, main pain points, and you 
end up with a smaller application that is easy to build, easy to verify, and low stakes if you completely break it! 

How does one create a mini-app? First, ask the oceanographers what are the pain points! They'll quickly tell you about the 
"dynamical core", which is where the basic physics are solved, i.e. the fluid mechanics part of the code. The crucial bit 
here is to adopt a "how hard can it be?" mentality, be positive. I was quickly pointed to the file `MOM_dyn_split_RK2.F90`, 
which contains the main RK2 timestepping procedure. Fortran is our friend here, there is not a lot of object orientation 
shenanigans happening, i.e. what you read is what you get (most of the time). The RK2 step can be summarized as:

```fortran 
!! The split RK2 scheme:
!!   1. Predictor phase:
!!      - Compute horizontal viscosity
!!      - Compute Coriolis/momentum advection (CAu, CAv)
!!      - Apply vertical viscosity
!!      - Advance barotropic mode (fast 2D dynamics)
!!      - Update layer thicknesses via continuity
!!   2. Corrector phase:
!!      - Recompute tendencies with updated state
!!      - Apply vertical viscosity
!!      - Final barotropic step
!!      - Final continuity update
!!
```

Don't get too distracted by initialization procedures, those might look ugly but once you pay attention they are simply setting arrays
to something. The above algorithm has done most of the work for us now! We now know what routines we need to look for:

- Horizontal/Vertical viscosities
- Coriolis 
- Barotropic solver 
- Continuity

A key concept here can be extracted from MOM, the first M means "modular". Each physics engine should be a module, i.e. it should be able to 
be run as standalone. Why? Because if I only want to optimize the Coriolis routines, why should I have to run everything else? So when creating a mini-app
always try to have a separation of concerns. 

We go into the magical place that is `MOM6/src/core`, which is where most of the routines we're interested in are. The viscosities are _parametrizations_, i.e. 
approximations to extremely annoying physics. If you're feeling adventurous code up a non-hydrostatic solver that tries to resolve the vertical physics. It is
_ a w f u l_. Here we have:

```
MOM_barotropic.F90
MOM_continuity_PPM.F90
MOM_CoriolisAdv.F90
```

So then we follow the `MOM_dyn_split_RK2.F90` to find the function names that are the most interesting, for simplicity we will do the `call btstep(..)`. Inside
a MOM6 subroutine, we can think of them like:

```fortran 
subroutine btstep(....)
real :: x
integer :: y
logical :: a 

! indices mapping 
is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke 

! bookkeeping (timing, profiling) 

if(OBC) then 

endif 

! halo updates (look for do_group_pass)
call create_group_pass(...)

! begin real work 
...

end subroutine btstep
```
This is not the rule everywhere, especially in helper functions and subroutines. But realizing that _most_ do some bookkeeping, accounting at the start, and then 
they check for configs, such as `OBCs` etc. will make your life less overwhelming. Once in, the magic is just realizing the loops are just loops! 

```fortran 
      do j=js,je ; do I=is-1,ie
        uhbt0(I,j) = uhbt(I,j) - find_uhbt(dt*ubt(I,j), BTCL_u(I,j)) * Idt
      enddo ; enddo
```

As you go into the routine you'll realize that there are patterns. Certain routines loop over certain indices, there's `u`s and `v`s and you start quickly 
noticing how the algorithms are structured. This is the way, you simply just continue doing this _with a lot of patience_. 

The mini-app is a work of hackathon and obsession, trying to compare lots of stuff (serial, openmp, openacc, cuda-fortran) just to get knowledge about the 
advantages and disadvantages of each approach. Testing performance without having to worry too much about breaking MOM6. 




## Navier Stokes -> stacked shallow water (adiabatic)
## Generalised vertical coordinates
## Vertical Lagrangian remapping
## Pressure forces
## Coriolis term
## Pressure solver - barotropic / baroclinic split
## Timestepping / Advection schemes
## Vertical velocity diagnostics
## MEKE
## EPBL
## KPP
## Submesoscale parameterisation
## Ice shelf code
## Ocean fresh-water forcing / water balance
## How to create a new MOM6 diagnostic
## Horizontal viscosity
## Advection schemes
## Bottom drag (Callum? Luwei?)
## Opacity for shortwave penetration
## For low res: GM / Redi eddy params
## Explicit tides, self attraction and loading
## Regional boundary conditions
## Bulk formulae
## Fortran 101 - control structure, where parameters are defined, some keywords e.g. private, intent, submodule, use module etc
## How to work out how a diagnostic was calculated (general syntax + keywords to search for)
## Different initialisation options - thickness + tracers, topography (all these preset options)
## Coord config vs use_regridding (layer vs ale mode) - how to change your vertical coordinate or target coordinate and what it means
## Regridding of diagnostics and vanished layer output - how to change in diag_table, what it means for output
## Equations of state and different types of salinity, temperature in initialisation, model running and output
## What does the sea ice model do for the ocean ?
## MOM5/MOM6 code differences?
## Shear mixing
## Internal tide mixing
## Lee wave mixing
## Bottom boundary mixing
## Background mixing
