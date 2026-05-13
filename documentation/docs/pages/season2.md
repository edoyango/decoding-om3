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


## Fortran 101 - control structure, where parameters are defined, some keywords e.g. private, intent, submodule, use module etc
Presenter: @edoyango (Edward Yang).

Date: 14/05/2026

Contributing to [MOM6](https://github.com/acCESS-nri/mom6) can be extra daunting if you're not used to programming in Fortran. These
notes aim to introduce Fortran to someone who might already be familiar with Python. And thankfully, most of the Fortran features
exercised in MOM6 have a Python equivalent. Here, we won't be looking at MOM6 code directly, because the code itself is quite long
-- even if the language features used aren't too complicated. Instead, we'll build a simple example program that uses some of the
concepts that MOM6 is built with. The program takes a 3d array (first two indices represent the lateral domain, and the last
represents the columns), and performs either a sum or max operation along the column, reducing the result to a 2d array that
represents the lateral domain only.

<details><summary>Python equivalent of the example program to be built</summary>

```python
import numpy as np


class Grid:
    def __init__(self, is_: int, js: int, ie: int, je: int, nz: int):
        self.is_ = is_
        self.js = js
        self.ie = ie
        self.je = je
        self.nz = nz


class ControlStructure:
    def __init__(self, initialized: bool = False, which_op: int = 1):
        self.initialized = initialized
        self.which_op = which_op


def do_something_along_column(cs: ControlStructure, g: Grid, arr1: np.ndarray) -> np.ndarray:
    if not cs.initialized:
        raise RuntimeError("Control structure not initialized!")

    if cs.which_op == 1:
        return _sum_along_column(g, arr1)
    if cs.which_op == 2:
        return _max_along_column(g, arr1)
    raise ValueError("Invalid operation provided! must be either 1 or 2")


def _sum_along_column(g: Grid, arr1: np.ndarray) -> np.ndarray:
    arr2 = arr1[:, :, 0].copy()
    for j in range(g.je - g.js + 1):
        for k in range(1, g.nz):
            for i in range(g.ie - g.is_ + 1):
                arr2[i, j] += arr1[i, j, k]
    return arr2


def _max_along_column(g: Grid, arr1: np.ndarray) -> np.ndarray:
    arr2 = arr1[:, :, 0].copy()
    for j in range(g.je - g.js + 1):
        for k in range(1, g.nz):
            for i in range(g.ie - g.is_ + 1):
                arr2[i, j] = max(arr2[i, j], arr1[i, j, k])
    return arr2

import numpy as np

from my_module import ControlStructure, Grid, do_something_along_column


if __name__ == "__main__":
    g = Grid(is_=1, js=2, ie=3, je=4, nz=5)
    cs = ControlStructure(initialized=True, which_op=1)

    input_array = np.ones((g.ie - g.is_ + 1, g.je - g.js + 1, g.nz))
    output_array = do_something_along_column(cs, g, input_array)

    print(np.sum(output_array))
```

</details>


### Programs and modules

Most of the code in MOM6 is organised into "modules" which usually relate to a certain area of the ocean physics. For example,
`MOM_barotropic.F90` contains the [barotropic solver code](http://github.com/ACCESS-NRI/MOM6/blob/2026.01/src/core/MOM_barotropic.F90),
and `MOM_tracer_advect.F90` contains the [tracer advection code](http://github.com/ACCESS-NRI/MOM6/blob/2026.01/src/tracer/MOM_tracer_advect.F90)
and so on. Modules contains code that can be reused in other modules or "programs". Modules cannot be directly compiled and run, and
so modules' code must be "used" from a program. The skeleton of this arrangement can look like:

```fortran
module my_module ! my_module is the name of the module
  implicit none
end module my_module

program my_program ! my_program is the name of the program
  use my_module
  implicit none
end program my_program
```

The module and program's boundaries in the code are deliniated by `program/module` and `end program/module` couples. The
program/module's name must follow the first `program/module` and matching `end program/module` (including the program/module name is
optional, but it's a requirement to include the name in modules in).

`implicit none` is a "quirk" of Fortran. It says that all variables' type must be declared. Otherwise, the compiler can make an
"educated guess" as to what the type it is, which can result in unexpected behaviour. Hence, it is best practice to always include
`implicit none` in every module and program.

### Subroutines and declaring variables

Programs can have runnable code. But as your codebase gets larger, it's likely that you will 1. want to organise the code in some
way to make it easier to understand and maintain (e.g. group code related to certain physics together), and 2. store code that is
reused in multiple places. Subroutines help facilitate this. Subroutines are similar to Python functions except that subroutines
don't return anything, and instead, they modify its arguments instead. Fortran does have functions also, but are used less
frequently in MOM6.

Like programs and modules, subroutines are bounded by `subroutine <name>` and `end subroutine <name>`. The subroutine's arguments
follow the name, followed by the type declaration of the arguments and local variables. Unlike Python, the types of all variables
must be declared and cannot change. In the example below, the variable attributes commonly used in MOM6 are shown (type, dimension,
and intent). Variables can be declared in any order, but a common convention (that MOM6 follows) is to declare the arguments first,
followed by local variables.

The example above is trivial as both the program or module have no code. So let's create our first subroutine (Fortran comments are
prefixed with `!`):

```fortran
module my_module

  implicit none

! says that subroutines/functions are declared after
contains

  ! this is declaring the subroutine's signature
  subroutine sum_along_column(is, ie, js, je, nz, arr1, arr2)

    ! each argument's type must be declared. Here we have:
    ! * type (integer/real)
    !   * real is equivalent to np.float32. However, MOM6 opts to control the precision at compile
    !     time.
    ! * dimension aka shape. No dimension means that variable is scalar. dimension(...) means the
    !   variable is an array.
    !   * dimension(a:b) means that for the given index, only indices a to b are defined.
    !     This a unique feature of Fortran that is utilised heavily in MOM6.
    ! * intent
    !   * `in`: the variable will only be read
    !   * `out`: the variable will be written to
    !   * `inout`: not shown here, but means that the variable may be read and/or written to.
    ! once an argument has been declared, it can also be used to declare others
    integer,                           intent(in)  :: is, ie, js, je, nz
    real, dimension(is:ie, js:je, nz), intent(in)  :: arr1
    real, dimension(is:ie, js:je),     intent(out) :: arr2
    ! all local variables must also be declared
    integer :: i, j, k

    do j=js,je
      ! copies can be done using array slicing
      arr2(:,j) = arr1(:,j,1) ! initialize the sum along columns
      ! MOM6 keeps nested loops on a single line, dilineated by colons.
      do k=2,nz ; do i=is,ie
        arr2(i,j) = arr2(i,j) + arr1(i,j,k) ! do the sum along columns
      enddo ; enddo
    enddo

  end subroutine sum_along_column

end module my_module
```

Note that subroutines are invoked with `call <subroutine>(...)`

### Array copies

Like with Python NumPy, arrays' contents can be copied between each other. If the arrays are of the same shape, they can be copied
with specifying array indices (`a = b`), or you may specify which slices to copy e.g. `a(1:10) = b(1:10)`, or `a(:, 1) = b(:)` etc.
Noting that Fortran accesses array elements/slices using round brackets `()` instead of square brackets `[]` common in other
languages. It is also worth noting that array assignments/copies copies the array's values (different from Python where `a = b`
means something different from `a[:] = b[:]`).

### Loops

Loops are dilineated by `do variable=start,end,step` and `enddo` - which is like `for variable in range(start,end+1,step):` in
Python. The main difference between fortran loop ranges and Python ranges is that the `end` is included in the range.

You may have also noticed that the loop ordering is a bit strange - where the outer loop iterates of the middle index, with the
outer index in the middle, and the inner-most loop iterating over the first. This is quite common in MOM6!

### Derived types

Derived types are similar to Python classes in that they can be instantiated and group related information. Like Python classes,
Fortran derived types has a name and attributes (or members). One of the key derived types in MOM6 is the
[`ocean_grid_type`](http://github.com/ACCESS-NRI/MOM6/blob/2026.01/src/core/MOM_grid_type.F90) which describe the grid extents
(including the computational and halo extents). It also stores other grid information like lateral dimensions of the columns,
masking etc. We can create a simple version of the grid type and use it in our subroutine:

```fortran
module my_module

  implicit none

  ! grid type to hold grid bounds                     Python analogue:
  type grid_type                                      ! class gridType:
    integer :: is, js   ! starting indices            !     def __init__(self, is, js, ie, je, nz):
    integer :: ie, je   ! ending indices              !         self.is = is ; self.js = js
    integer :: nz       ! number of vertical layers   !         self.ie = ie ; self.je = je
  end type grid_type                                  !         self.nz = nz

contains

  ! we can replace our grid indices with a grid_type
  subroutine sum_along_column(g, arr1, arr2)
    ! We can still use the grid_type's members to define subsequent variables
    ! below, `g` is an instance of the `grid_type`. This `g`'s members are being accessed with `%`,
    ! which in python would be `g.is`. `g`'s members are then used to size the input/output arrays.
    type(grid_type),                             intent(in)  :: g    ! g is an instance of grid_type
    real, dimension(g%is:g%ie, g%js:g%je, g%nz), intent(in)  :: arr1
    real, dimension(g%is:g%ie, g%js:g%je),       intent(out) :: arr2
    integer :: i, j, k

    do j=g%js,g%je
      arr2(:,j) = arr1(:,j,1) ! initialize the sum along columns
      do k=2,g%nz ; do i=g%is,g%ie
        arr2(i,j) = arr2(i,j) + arr1(i,j,k) ! do the sum along columns
      enddo ; enddo
    enddo

  end subroutine sum_along_column

end module my_module
```

Let's also introduce another pattern used in MOM6: the "control structure". Each module will have its own control structure that
mostly stores information to control MOM6's behaviour e.g. which algorithm to use or whether a certain physics is turned on or not.
In our simple example, the control structure will simply control whether to do a sum or max along columns.

```fortran
module my_module

  implicit none

  private ! this says that by default, contents of this module aren't accessible.

  !< grid type to hold grid bounds
  type grid_type
    integer :: is, js   !< starting indices
    integer :: ie, je   !< ending indices
    integer :: nz       !< number of vertical layers
  end type grid_type

  !< control structure
  type control_structure_type
    logical :: initialized = .false. !< whether the control structure is initialized - defaults to .false.
    integer :: which_op = 1          !< which operation to do - default is 1 (sum)
  end type control_structure_type

  ! Explicitly say which things can be used - do_sum/max_along_column cannot be directly used
  public :: grid_type, control_structure_type, do_something_along_column

contains

  !< this subroutine either does a max or sum along columns
  subroutine do_something_along_column(cs, g, arr1, arr2)
    type(control_structure_type),                intent(in)  :: cs   !< control structure
    type(grid_type),                             intent(in)  :: g    !< grid type
    real, dimension(g%is:g%ie, g%js:g%je, g%nz), intent(in)  :: arr1 !< input array
    real, dimension(g%is:g%ie, g%js:g%je),       intent(out) :: arr2 !< output array

    if (.not.CS%initialized) error stop "Control structure not initialized!"

    if (CS%which_op == 1) then
      call sum_along_column(g, arr1, arr2)
    elseif (CS%which_op == 2) then
      call max_along_column(g, arr1, arr2)
    else
      error stop "Invalid operation provided! must be either 1 or 2"
    endif

  end subroutine do_something_along_column

  !< performs a sum reduction along columns
  subroutine sum_along_column(g, arr1, arr2)
    type(grid_type),                             intent(in)  :: g    !< grid type
    real, dimension(g%is:g%ie, g%js:g%je, g%nz), intent(in)  :: arr1 !< input array
    real, dimension(g%is:g%ie, g%js:g%je),       intent(out) :: arr2 !< output array
    integer :: i, j, k ! loop indices

    do j=g%js,g%je
      arr2(:,j) = arr1(:,j,1) ! initialize the sum along columns
      do k=2,g%nz ; do i=g%is,g%ie
        arr2(i,j) = arr2(i,j) + arr1(i,j,k) ! do the sum along columns
      enddo ; enddo
    enddo

  end subroutine sum_along_column

  !< performs a max reduction along columns
  subroutine max_along_column(g, arr1, arr2)
    type(grid_type),                             intent(in)  :: g    !< grid type
    real, dimension(g%is:g%ie, g%js:g%je, g%nz), intent(in)  :: arr1 !< input array
    real, dimension(g%is:g%ie, g%js:g%je),       intent(out) :: arr2 !< output array
    integer :: i, j, k ! loop indices

    do j=g%js,g%je
      arr2(:,j) = arr1(:,j,1) ! initialize the sum along columns
      do k=2,g%nz ; do i=g%is,g%ie
        arr2(i,j) = max(arr2(i,j), arr1(i,j,k)) ! do the max along columns
      enddo ; enddo
    enddo

  end subroutine max_along_column

end module my_module
```

The module is much larger now - the control structure type has been added and two subroutines have also been added.
`max_along_column` is almost identical to `sum_along_column` except that it does a `max` operation instead of `+`.

### If statements

In the module, we've also added a 3-branch if statement. If statements look quite similar to Python's except that the evaluation
must be put it into brackets and is followe by `then`: i.e. `if (statement) then`. Otherwise the semantics are identical.

### Module public and private

To hide details, MOM6 likes to leverage "public" and "private" statements in modules. MOM6 modules will declare everything as
`private` by default (by having an unqualified `private` clause), and then explicitly list the objects that should be accessible
outside of the module with `public :: list, of, objects, and, procedures`. Private things are visible to other things within the
module, but not outside.

### Documenting comments

MOM6 uses [Doxygen comments](https://www.doxygen.nl/manual/docblocks.html) that automatically generate documentation for the code.
These type of comments are sentineled with !< (as opposed to only !). Procedures (functions and subroutines), types and members, and
arguments must be documented.

### Using modules

Let's finish our program and use the module code! The program will be very simple -- uses hardcoded values to initialize everything.

```fortran
program my_program
  ! specify which things we want from the module
  ! Python analogue: from my_module import grid_type, control_structure_type, do_something_along_column
  ! note that because of the privacy in the module, `sum_along_column` is not included in the
  ! `public :: ...` list. This means that importing or calling `sum_along_column` from this
  ! program would fail. Adding `sum_along_column` to the `public :: ...` list would make it
  ! available here.
  ! Note that excluding `only` would lead to everything public being available. However, MOM6
  ! ensures to make imports explicit. This is, in general, good practice as it makes it easier to
  ! identify where something used inside the program/module comes from.
  use my_module, only: grid_type, control_structure_type, do_something_along_column
  implicit none
  ! variable declaration
  type(grid_type) :: g
  type(control_structure_type) :: cs
  ! use allocatable arrays for dynamic array sizes in the program
  real, allocatable :: input_array(:, :, :), output_array(:, :)

  ! derived type instances don't have to be constructed (the construction is implied in the declaration)
  ! but they can be initialized with a default type constructor, or simply setting the members.
  g = grid_type(is=1, js=2, ie=3, je=4, nz=5)
  cs%initialized = .true.
  cs%which_op = 1

  allocate( &
    input_array(g%is:g%ie, g%js:g%je, g%nz), &
    output_array(g%is:g%ie, g%is:g%ie) &
  )

  input_array = 1.0

  call do_something_along_column(cs, g, input_array, output_array)

  ! print the sum of the array to the terminal
  write(*, *) sum(output_array)

end program my_program
```

A program itself looks quite similar the subroutines above, except there aren't any arguments. Near the top we `use` the module and
select which things we want. This convention isn't mandatory, but is closely followed by MOM6.

### Allocatable arrays

Unlike subroutines, programs cannot get array sizes by simply passing them in. So, often arrays are made given the "allocatable"
attribute, which lets the program "allocate" the array based on some user input or similar. Here, we allocate the arrays' based on
hardcoded values, but that could be changed to use CLI arguments or something.

### Compiling and running the example

To compile the program, you can put the above module and program into the same file, say `example.f90`, and compile it with
`gfortran example.f90 -o example.x`. You can then execute it with `./example.x`. You should get `45` with many zeros printed!

MOM6 is a much more complex codebase with many more dependencies and consequently more complex build system.


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
