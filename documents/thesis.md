**A FDTD Model of a 2D Kirchhoff Thin Plate Model and Coupled Strings**

![image](../FigsGlobal/eushield_Colour.pdf){width="4.5cm"}

A final project dissertation submitted in partial fulfilment\
of the requirements for the degree of\

**Master of Science (MSc)**\

***Acoustics and Music Technology***\

Acoustics and Audio Group\

Edinburgh College of Art\
University *of* Edinburgh\

August 2017

Supervisor 1: Dr Stefan Bilbao\

<span> </span>

<span> </span>

Declaration {#declaration .unnumbered}
===========

I do hereby declare that this dissertation was composed by myself and
that the work described within is my own, except where explicitly stated
otherwise.

Matthew Hamilton\
August 2017

<span> </span>

Acknowledgements {#acknowledgements .unnumbered}
================

For my wife Helen, without whose help and support I would not have
achieved any of the following.

<span> </span>

<span> </span>

<span> </span>

Introduction {#chapter1}
============

With ever-increasing computer processing power and its availability to a
home user, the feasibility of implementing more involved and
processor-intensive computation increases. With respect to computer
music, the potential for realising new and interesting ideas is far
greater than it was even 10 years ago. This means the ability to create
much more complex and intricate ways of processing and interacting with
sound. One such way of generating sound is through physical modelling
synthesis. In the past, physical modelling has been less feasible as the
methods involved tended to be computationally intensive. Now these same
methods are not only viable but they are ever more capable of being
achieved in real-time.

Physical-modelling in this sense is the description of the movement of
some physical system to a digital domain. The equations of motion which
are used to predict the behaviour of musical systems, such as membranes
and plates strings &c.. are used to generate sound. The benefits are
that the sounds tend to be richer than other modes of sound synthesis as
well as the fact that they are more closely informed by the behaviour of
real instruments. One of the ways of physical modelling of musical
systems is through use of finite difference time domain method (FDTD).
FDTD is a method of translating the differential operators of partial
difference equations (PDEs) used in equations of motion. As such, they
are a little more tricky to get to grips with than traditional methods
of sound synthesis. Should one wish to begin creating FDTD systems to
make instruments or audio effects it can help to understand some of the
underlying theory involved in FDTD.

The purpose of this paper is to explore the means for creating more
complex FDTD instruments by coupling, in particular coupling a plate
model with string models. This will allow the creation of
chordophone-like instruments that, depending on the configuration,
resemble either a piano, a harp or a guitar.

This chapter will deal with basics of FTDT, the operations involved and
some methods for beginning to tackle this kind of problem. The second
chapter will tackle the thin plate and and its respective FTDT model. It
will also discuss making the model more complex, what this involves and
how to go about deriving the correct FDTD scheme. Chapter 3 will concern
coupling the FDTD plate model with, firstly, a 1D wave and then a stiff
string model. Finally Chapter 4 concerns the additional measure that can
be taken to produce more interesting sounds as well other ways of
interacting with the model and extensions that can be pursued beyond
what is covered in this paper.

If unacquainted with the type of algebra contained within this paper,
explicit derivations have been provided in Appendix A and are referenced
throughout. Appendix B provides some example function and FDTD schemes
which should aid in simplifying the coding process.

FDTD Theory
-----------

First, it is worthwhile giving an overview of FDTD operators and how
they work. FDTD operators represent an approximation of differential
operators in continuous time and space. This centres around a difference
between two points in space or time divided by the meaasured distance
between. Since FDTD is in the time domain the goal is to find an unknown
displacement that is forward in time. $$\begin{aligned}[c]
\delta_{x+} & \triangleq \frac{1}{h} ( e_{x+} - 1 )\\
\delta_{x-} & \triangleq \frac{1}{h} ( 1 - e_{x-} ) \\
\delta_{x\cdot} & \triangleq \frac{1}{2h} ( e_{x+} - e_{x-} )
\end{aligned}
\qquad
\begin{aligned}[c]
\delta_{t+} & \triangleq \frac{1}{k} ( e_{t+} -1  ) \\
\delta_{t-} & \triangleq \frac{1}{k} ( 1 - e_{t-} ) \\
\delta_{t\cdot} & \triangleq \frac{1}{2k}( e_{t+} - e_{t-} )
\end{aligned}$$

Operators behave the same way in space as they do in time and involve
the difference between two points which are shifted and divided by the
spacing that separates them. In terms of differential term
$\frac{dy}{dt}$ the $\frac{1}{dt}$ can be considered the continuous form
of 1/k. For the time operators, $\delta_{t}$, the variable k is the
period of time between samples (i.e. 1/sampling rate). The smaller k is,
the more accurate a representation it can be considered of the
continuous case. Operators can be combined to create approximations of
higher order difference. $$\begin{aligned}
\delta_{xx} &= \frac{1}{h} ( \delta_{x+} - \delta_{x-} ) = \frac{2}{h} ( \delta_{x\cdot} - \delta_{x-} )\\
\delta_{xxxx} &= (\delta_{xx})(\delta_{xx})\end{aligned}$$

It is hopefully apparent how each operator functions; if not,
[@BilbaoNSS] is an excellent source, especially with its perspective of
FDTD methods for musical applications. More complex operators such as
the laplacian, $\Delta$ and the biharmonic $\Delta\Delta$ will be
discussed later. It will hopefully show that spatial difference can be
mostly automated and temporal difference will be restricted to
concerning three time steps.

Notation
--------

The notation for this paper will be in terms of dot and primes for
temporal and spatial derivatives respectively. Given that the plate is
two-dimenional, there is little overlap in the notation between systems.
This should hopefully assist in minimising any confusion.

For a 2D FDTD scheme, of greater concern are the laplacian and
biharmonic operators.

$$\begin{aligned}
\Delta = \left[\pdv[2]{}{x} + \pdv[2]{}{y} \right], \quad \Delta\Delta = \left[\pdv[4]{}{x} + 2\frac{\partial^{4}}{{\partial^{2}x}{\partial^{2}y}} +  \pdv[4]{}{y} \right] = \left[\pdv[2]{}{x} + \pdv[2]{}{y} \right]^2
\label{eq:2Dops}\end{aligned}$$

The laplacian is a partial second order term in each of the spatial
dimensions and the biharmonic is simply the laplacian squared. The
implications this has on a FDTD scheme will be discussed later. It is
worth familiarising oneself with the terms in equation \[eq:2Dops\]. If
new to these terms, the notation can do more to obfuscate the meaning
than make it clearer. In other texts $\nabla^2$ may be used instead of
$\Delta$, the latter will be used in this paper.

The FDTD systems that will be considered are the plate and the string.
These will always be represented by ‘u’ for the plate and ‘w’ for the
string. These terms both represent a vector of displacements/amplitudes.
The finer details will be discussed later. The superscript for these
terms refers to time step, ‘n’ and the subscript the spatial index. For
instance, $u_{n+1}$ refers to the state of the system one sample ahead
of the current time step. Given the approach that will be taken to the
discussed schemes, the spatial index will rarely be used as systems for
the plate and string will contained in the vector form.

FDTD Matrices {#chapter1:sec_LinDex}
-------------

![The transformation of a 2D grid to a linear indexed vector. The grid
co-ordinates are aligned with their corresponding linear index value.
$N_y$ is the number of grid points in the y-axis.<span
data-label="figs:LinDex"></span>](../Chapter_1/_Figs/LinDex){width="8cm"}

As had been previously stated, the end goal for a FDTD is to compute an
unknown value, forward in time, in terms of previously computed values.
In the confines of this paper there are three time steps that are
considered. The forward step that needs computed, $u^{n+1}$, the current
time step $u^n$ and previous step $u^{n-1}$. The FDTD will be organised
in vectors with linear indexing. Even though it might seem the obvious
first choice for a 2D model to be represented as a grid, a vector is
more readily operated on by matrices. The organisation of linear
indexing for the plate is given in Figure \[figs:LinDex\]. The spatial
differences can be contained within $N \times N$ matrices, where N is
the number of grid points in the system. The update at each time step to
find the next state will always be in the form

$$Au^{n+1} = Bu^n + Cu^{n-1}
\label{eq:FDupdate}$$

In equation \[eq:FDupdate\], A, B and C refer to coefficient matrices
consiting of the spatial FDTD operations that relate to each time step.
Building these matrices will be the main task for creating the schemes
discussed here. The benefits of sticking to this form are that very
little change should actually be required when altering and adding
greater complexity to the schemes, which will hopefully be conveyed
throughout the remainder of the paper.

The same thinking can be applied to spatial difference operators. As an
example, take the definition of the forward difference operator for the
1D wave $\delta_{\eta+}$.

$$\delta_{\eta+} w = \frac{1}{h^{2}} ( u_{\eta+1} - u_{\eta} )$$

Carrying out this operation across the whole domain of w can be
represented as the following matrix operation:

$$\frac{1}{h^{2}}D_{\eta+}\begin{bmatrix}
               w_{1} \\[0.3em]
               \vdots \\[0.3em]
               w_{N-1}\\[0.3em]
               w_{N}
             \end{bmatrix} 
 = \frac{1}{h^{2}} \begin{bmatrix}
             -1 & 1 & &   \\[0.3em]
               & \ddots & \ddots &    \\[0.3em]
                & & -1 & 1 \\[0.3em]
                & &  & -1
             \end{bmatrix} \begin{bmatrix}
               w_{1} \\[0.3em]
               \vdots \\[0.3em]
               w_{N-1}\\[0.3em]
               w_{N}
             \end{bmatrix}$$

Building matrices in this manner simplifies the problem of higher order
FD spatial operators. As another example for 1D wave,
$\delta_{\eta\eta} = \frac{1}{h}(\delta_{\eta+} - \delta_{\eta-}) $

$$\begin{aligned}
\delta_{\eta\eta} & = \frac{1}{h^2}(D_{\eta+} - D_{\eta-}) \\
 & = \frac{1}{h^{2}} \begin{bmatrix}
             -1 & 1 & &   \\[0.3em]
               & \ddots & \ddots &    \\[0.3em]
                & & -1 & 1 \\[0.3em]
                & &  & -1
             \end{bmatrix}
		\begin{bmatrix}
               1 &  & &   \\[0.3em]
              -1 & 1 & & \\[0.3em]
                  & \ddots & \ddots &    \\[0.3em]
                  & & -1 & 1
             \end{bmatrix} \\ 
             & = \frac{1}{h^{2}} \begin{bmatrix}
               -2 & 1 & &   \\[0.3em]
               1 & -2 & 1 & \\[0.3em]
                  & \ddots & \ddots & \ddots  \\[0.3em]
                  & & 1 & -2
             \end{bmatrix} = \frac{1}{h^{2}} D_{\eta\eta}\end{aligned}$$

This approach turns the construction into a trivial problem that
requires only a simple matrix multiplied by the correct coefficient.
Approaching the problem in this manner is a much more attractive method
than the linear equivalent. This modular approach allows for much
greater separation of the various elements in the scheme in comparison
to coding it explicitly. The only extra consideration is the
modifications to the matrix coefficients as a result of boundary
conditions. This will be covered in section \[chapter2:sec\_boundaries\]
though the functions in \[appB:sec1\] provide one approach to automating
the solution. It is recommended to create your own FDTD matrix
functions; this will be a valuable aid in streamlining any code and will
provide a handy tool that can be reused. Though API and toolkits are
readily available, coding by hand should bring to light any
uncertainties or misconceptions. Functions also help to reduce any
needless repetition and improve code readability.

Kirschoff Thin Plate {#chapter2}
====================

The first stage of coupling two FD schemes together is having something
to couple to. It is best to begin with the thin plate model as it is the
more tricky of the two schemes. It should also help to clarify the
processes involved in making an FDTD scheme in general. The plate model
that will be considered is the Kirschhoff plate model. The model and
accompanying finite difference scheme will be stated strictly in
cartesian coordinates. It is of course possible to create a FDTD scheme
in radial co-ordinates and in fact, this is probably the first step to
creating a circular membrane or shell model such as for drums and
cymbals. For more on radial coordinates and plates see [@BilbaoNSS
Section 12.6]. The derivation of the these equations of motion will not
be discussed in any great detail in this paper. This has been covered
succinctly elsewhere, [@rossing2004principles] has been very helpful for
getting to grips with the physical analysis. The FDTD model of the plate
that will be dealt with here is that in [@BilbaoNSS] and by extension
the physical model from [@morse1968theoretical].

Thin Plate Equation {#chapter2:sec1}
-------------------

The Kirschhoff plate model, equation \[eq:kirschhoffplate\], is in
essence an extension of the 1D bar equation [@rossing2004principles].

\[chapter2:sec1:eq1\]
$$\rho H \ddot{u}  = -D\Delta\Delta u, \quad D = \frac{E H^3}{12(1- \nu)} 
\label{eq:kirschhoffplate}$$

On the left hand side the displacement of the plate in this case is $u$,
denisty $\rho$, plate thickness $H$.On the right is the fourth order
term, the biharmonic operating on u and the constant D containing the
Young’s Modulus of the material $E$ and the poisson ratio $\nu$.
Bringing over the density and thickness gives us a more simplistic
equation \[eq:kirschhoffplatesimple\] with one constant $\kappa$ which
will be easier to implement in a FDTD scheme. It is advisable, when
adding extra terms, to begin afresh from \[eq:kirschhoffplate\],
especially when including a force term. It is easy to forget the
inclusion of the $\rho H$ term and the scheme will be incorrect without
it.

$$\ddot{u}  = -\kappa^{2}\Delta\Delta u, \quad \kappa = \sqrt{\frac{E H^2}{12\rho(1- \nu)} }
\label{eq:kirschhoffplatesimple}$$

In this model the plate is isotropic, which is to say that wave velocity
is equal in all directions. Anisotropy is encountered in materials like
wood and will be discussed in Chapter \[chapter4\].

Thin Plate Finite Difference Model {#chapter2:sec2}
----------------------------------

A FDTD model for thin plate is given in \[eq:kirschoffFDTD\]. It is a
direct translation of the operators in \[eq:kirschhoffplate\].

$$\rho H \delta_{tt}u  = -D\delta_{\Delta\Delta} w \quad \to \quad \delta_{tt}u  = -\kappa^{2}\delta_{\Delta\Delta} w
\label{eq:kirschoffFDTD}$$

Deriving an update that can be used in a scheme is shown in equation
\[eq:kirschoffFDTDupdate\].

$$\begin{aligned}
\rho H \delta_{tt}u &= -D\delta_{\Delta\Delta}u \nonumber \\
\delta_{tt}u &= -\kappa^{2}\delta_{\Delta\Delta}u\nonumber \\
\frac{1}{k^2}(u^{n+1} - 2u^{n} + u^{n-1}) &= -\kappa^{2}\delta_{\Delta\Delta}u\nonumber \\
u^{n+1} &= -k^2\kappa^{2}\delta_{\Delta\Delta}u + 2u^{n} - u^{n-1} \nonumber \\
&= -\frac{k^2\kappa^{2}}{h^4}D_{\Delta\Delta}u + 2u^{n} - u^{n-1} \nonumber \\
 &= -\mu^2 D_{\Delta\Delta}u + 2u^{n} - u^{n-1} \nonumber \\
\label{eq:kirschoffFDTDupdate}\end{aligned}$$

The addition here is the scheme variable $\mu^2$ which is defined as:-

$$\mu = \frac{k\kappa}{h^2}$$

Here, $k$ is the time step and $h$ is the grid spacing which is to the
power 4, given that $\Delta\Delta$ consists of fourth order terms. To
set the plate in motion, some initial conditions, in the form of a
raised cosine, at some point on the plate can be used (see Appendix
\[cd:plateFD\] lines 98-102 and 118-119).

When beginning to solve for a FD update the benefits of the matrix
approach become more apparent. Take for instance the update for the the
biharmonic written out explicitly.

$$\begin{aligned}
\delta_{\Delta\Delta} & = \frac{1}{h^{4}} (e_{x+1} - 2 + e_{x-1} e_{y+1} - 2 + e_{y-1})^{2} \nonumber \\
 & =  \frac{1}{h^{4}}(20 - 8(e_{x+} + e_{x-} + e_{y+} + e_{y-})  + ... \nonumber \\
&\quad (e_{x+,y+} + e_{x+,y-} + e_{x-,y+} + e_{x-,y-}) + ... \nonumber \\
&\quad (e_{x+2} + e_{x-2} + e_{y+2} + e_{y-2}))\end{aligned}$$

This term would need to be applied to every point on the plate, not
having taken into account boundary conditions, which is a little
unwieldy. In matrix form the aesthetics are a little less distressing.
Moving forward in the form suggested in \[eq:FDupdate\],
\[eq:kirschoffFDTDupdate\] can be broken down into three matrix
operations.

$$\underbrace{(I)}_{A}u^{n+1}  = \underbrace{(-\mu{2}D_{\Delta\Delta} +2I}_{B})u^{n}  -\underbrace{(I)}_{C} u^{n-1}
\label{eq:PlateMatNoLoss}$$

The three matrices, A, B and C operate on the vector u, that is the thin
plate in this case. The matrices ‘A’ and ‘C’ may seem trivial in this
form but, as greater complexity is added these become correspondingly
more complex. It is worth noting at this point that the ‘$-$’ sign in
front of the C matrix can be included, $Au^{n+1} = Bu^{n} + Cu^{n-1}$ or
kept outside $Au^{n+1} = Bu^{n} - Cu^{n-1}$. For this paper the first
approach will be taken.

### Boundary Conditions {#chapter2:sec_boundaries}

This paper will primarily consider simply supported and clamped boundary
conditions. In particular, clamped conditions, as these should closely
relate to the behaviour of a soundboard of a string instrument as
suggested in [@GoughColin2015Vpm]. Since the plate equation is forth
order in both the x-axis and y-axis, there are two conditions in each
set. Both sets have a Dirichlet condition where the boundary points are
fixed at zero. The Clamped set has a Neumann condition that the gradient
is equal to zero. Simply Supported has a condition that the curvature is
fixed at 0. In FDTD this can be represented as:

$$\begin{aligned}
u_{0} & = u_{N}  =  0 \label{eq:con1}  \\ 
\delta_{x\cdot}u_{0} & = \delta_{x\cdot}u_{N} = 0 \label{eq:con2}  \\ 
\delta_{xx}u_{0} & = \delta_{xx}u_{N} = 0 \label{eq:con3}\end{aligned}$$

The implications of these conditions are not immediately apparent. In
fact, in second order equations \[eq:con2\] and \[eq:con3\] are not
required. In fourth order, however, we now have to take into account a
‘ghost’ point, $u_{-1}$, when considering locations neighbouring a
boundary. As an example, take the following forth order operation on a
point $u_{1}$ where the domain is between 0 and N

$$\delta_{xxxx} u_{0}  = \frac{1}{h^2}\left(u_{3} - 4u_{2} + 6u_{1} - 4u_{0} + u_{-1} \right)$$

The point $u_{-1}$ is non-existent, but, using the other condition we
can derive an expression in terms of known points.

$$\begin{aligned}
u_{-1} & = -u_{1} \quad \text{[Clamped]}  \label{eq:bound1} \\
u_{-1} & = u_{1} \quad \ \ \text{[Simply Supported]} \label{eq:bound2} \end{aligned}$$

For derivation see Appendix A, Equation \[dvn:boundaries\]

### Constructing the Biharmonic

As illustrated in equation \[eq:2Dops\], the biharmonic is simply the
laplacian operator squared. The same can be extended to the equivalent
FDTD operators. There are a number of ways to approach this problem.
Hard coding will impose greater limitations and not recommended
although, it can be a useful learning exercise if the matrices are is
kept small. The best means of constructing the biharmonic is
procedurally.

$$\begin{aligned}
\delta_{\Delta\Delta} = \delta_{\Delta}\delta_{\Delta} = (\delta_{xx} + \delta_{yy})^2\\
\frac{1}{h^4}D_{\Delta\Delta} = (\frac{1}{h^2}D_{\Delta})\frac{1}{h^2}D_{\Delta} = (\frac{1}{h^2}D_{xx} + \frac{1}{h^2}D_{yy})^2
\label{eq:2DFDop}\end{aligned}$$

Constructing the biharmonic, beginning from the $D_{xx}$ and $D_{yy}$
matrices, allows for the easy inclusion of boundary conditions on the
plate. For clamped and simply supported conditions, this is just a case
of altering two points in the second order matrices.

$$\begin{aligned}
D_{yy} & = \frac{1}{h^{2}} \begin{bmatrix}
               -2 & b & &   \\[0.3em]
               1 & -2 & 1 & \\[0.3em]
                  & \ddots & \ddots & \ddots  \\[0.3em]
                  & & b & -2
             \end{bmatrix}
             \label{eq:boundmat}\end{aligned}$$

In \[eq:boundmat\], changing the constant $b$ to $0$ will implement
simply supported conditions, and $-2$, clamped conditions. Both of these
alterations are determined by \[eq:bound1\] and \[eq:bound2\]. Since the
plate grid is indexed linearly, there can be technical and conceptual
hurdles to overcome in constructing the $D_{xx}$ matrix. In the MATLAB
environment tools, such as the `kron()` function, prove very useful. In
the case of the function in Appendix \[cd:biharm\] that creates the
biharmonic, this is the approach that has been taken. Creating these
matrices individually can also be useful for implementing any FDTD
scheme (see Appendix \[cd:fidimat\])

An excellent resource and further reading can be found in [@ATorin_PhD],
which has some very helpful illustrations for beginning to understand
the vectorisation of grids and the construction of the matrices for FDTD
plates.

![Laplacian ($D_\Delta$) and Biharmonic ($D_{\Delta\Delta}$) Matrix
Patterns<span
data-label="figs:SpyMats"></span>](../Chapter_2/_Figs/SpyMats){width="10cm"}

![Biharmonic Stencil With linear index offsets for x and y axis. P
varies with number and type of boundaries it is adjacent to. $N_y$
represents the number of grid points in the y-axis<span
data-label="figs:BiharmStencil"></span>](../Chapter_2/_Figs/BiharmStencil2){width="10cm"}

Figure \[figs:SpyMats\] shows the pattern that should be obtained from a
biharmonic and laplacian matrix. In the case of Figure \[figs:SpyMats\]
the boundary points that are fixed at 0 have been left. In the case of
sparse matrix this should not add any extra strain to computation, but
it would be advisable to remove these calculations in other
implementations. Using linear indexing, the outer diagonals are the
coefficients for differences in the x-axis, while the y-axis differences
are those points adjacent to the diagonals.

Figure \[figs:BiharmStencil\] illustrates the stencil for the biharmonic
matrix. The value of the central point P depends on the number and type
of boundaries that are adjacent to it. All other points either keep the
same value or are zero, depending on grid point index. Correctly zeroing
out these points is given in Appendix \[cd:biharm\] lines 29-30.

For linear indexing the correct index for each point needs to be
modified by the given values on the x and y axis. For simply supported
conditions P is 20 for internal grid points, 19 when adjacent to a
single boundary and 18 when in a corner. The pattern for these points in
clamped conditions is 20, 21 and 22 respectively. There is not too much
effort to be had in mixing boundary conditions but for the purposes of
this paper it has be restricted to being unanimous across all
boundaries. The alteration in \[eq:boundmat\] should hopefully show that
enforcing arbitrary conditions for each boundary is fairly trivial.

### Energy Analysis

![Motion of a thin plate FD model with frequency dependant loss under
clamped conditions and initial conditions of a raised cosine.)<span
data-label="figs:plateWave"></span>](../Chapter_2/_Figs/plateWave){width="12cm"}

Energy Analysis is a good means of sanity checking and proof that a FD
model is showing stable behaviour. This can be very help full when
fault-finding during coding, particularly with the addition of greater
complexity. For example, should a model be stable when lossless and
unstable when loss conditions are added, it is likely that the
instability stems from the implementation of loss conditions.

Deriving energy follows the form of multiplying by the velocity and
integrating over the domain. This will provide the kinetic and potential
energy terms along with the boundary conditions $\mathcal{B}$ that need
to be satisfied. In a lossless system, the change in energy
$\frac{dE}{dt}$ should equal 0. For the thin plate model, the change in
energy is defined in equation \[eq:Energy\]. The total energy $E$ has
been marked. The derivation for this can be found in Appendix A
\[eq:PlateEnergy\]

$$\label{eq:Energy}
\frac{\partial}{\partial t}\underbrace{\left[\frac{1}{2}\|\dot{u}\|^{2}_{\mathcal{D}} + \frac{\kappa^2}{2}\|\Delta u\|^{2}_{\mathcal{D}}\right]}_{E}  = 0$$

For the discrete case we sum over the domain as there are a limited
number of points on the grid. In this case, there is no longer an
infinitesimally small distance being integrated over but rather a
measurable quantity in the form of the grid spacing. The energy in the
discrete case is given by equation \[eq:DiscEnergy\] (also see Appendix
\[dvn:energyDisc\]). In this instance the grid spacing is for x and y
are the same (i.e. $h_x = h_y = h$) There is a slight alteration that
needs to made when deriving this in discrete time. The change stems from
summation by parts(see [@BilbaoNSS Section 5.2.10]).

$$\label{eq:DiscEnergy}
\delta_{t+}\left[\frac{h^2}{2k^2}(u^n-u^{n-1})^T(u^n-u^{n-1}) + \frac{\kappa^2}{2h^2}(D_{\Delta} u^n)^T(D_{\Delta} u^{n-1})\right] =0$$

When implementing this in a FDTD scheme, only a very short time period
need be considered, which should save computation. The result should be
something similar to that found in Figure \[figs:dEdt\], the scattering
pattern is created by rounding error as the difference in values is near
machine epsilon.

![Rounding Error Pattern in change in energy($\frac{dE}{dt}$)<span
data-label="figs:dEdt"></span>](../Chapter_2/_Figs/dEdt){width="10cm"}

### Unit Checking

As the FDTD schemes increase in complexity, especially when coupling two
systems, it is easy to get lost in a jumble of coefficients. For the
most part, it tends to get uglier before it gets any better. Unit
checking is always a good sanity check for any fault finding. No units
in these schemes are scaled and the units of each term should be equal.

$$\rho H \pdv[2]{}{t}\left[u\right]  = -D\Delta\Delta u, \quad D = \frac{E H^3}{12(1- \nu)}$$
Considering just the units on the LHS
$$\rho = kgm^{-3}, \quad H = m, \quad \pdv[2]{}{t} = s^{-2}, \quad w = m$$

All terms on the LHS combine into a single term in units of
$kgm^{-1}s^{-2}$. So long as the RHS reflects the same units then we
know that we can proceed without issue. For the term $D$, the poisson
ratio is without units. We only consider the Young’s modulus and the
height cubed. The displacement u in $m$ and the height $H$ in $m^{3}$
cancel out the biharmonic which is units of $m^{-4}$. All that is left
is $E$ which is units of $kgm^{-1}s^{-2}$ which matches the units on the
LHS. This is a very useful tool for sanity checking, especially when
force is introduced, as the equations can start to become very
cluttered.

Loss
----

So far, the system does not have very much musical function as it will
continue moving forever. By adding loss, a decay in the energy will be
introduced, creating a sound of limited duration. In the first instance
a generic loss term can be added to \[eq:kirschhoffplate\], gives
equation \[eq:pWgloss\].

$$\begin{aligned}
  \rho H \ddot{u} &= -D\Delta\Delta u -2\rho H \sigma_{0}\dot{u} \nonumber  \\
   \ddot{u} &= -\kappa^2\Delta\Delta u -2\sigma_{0}\dot{u}
   \label{eq:pWgloss}\end{aligned}$$

Below is a FDTD scheme and its corresponding update for \[eq:pWgloss\]
(see Appenix \[dvn:plateFDTDWloss\]).

$$\begin{aligned}
\label{eq:pFDTDWloss}
\rho H \delta_{tt}u &= -D\delta_{\Delta\Delta}u - 2\rho H \sigma_{0}\delta_{t\cdot}u \nonumber \\
\underbrace{(1+ k^{2}\sigma_{0})}_{A}u^{n+1} &= \underbrace{(-\mu^2 D_{\Delta\Delta} + 2I)}_{B}u^{n} - \underbrace{(1 - k^{2}\sigma_{0})}_{C}u^{n-1}\end{aligned}$$

The A and C matrices are no longer empty. These matrices will primarily
deal with loss terms in the following FDTD schemes. Since A is a
diagonal matrix, its values can be precomputed and integrated into the B
and C matrices.

Applying Force {#chapter2:sec_force}
--------------

Rather than applying an initial displacement, a force can be added to
the plate allowing for multiple excitations. This allows for interaction
with the model and the ability to start using it like a musical
instrument.

$$\label{plateWForce}
\rho H \ddot{u} = -D\Delta\Delta u -2\rho H \sigma_{0}\dot{u} + \delta(x-x_0, y-y_0)f(t)$$

In \[plateWForce\] the force is applied with a dirac delta. The dirac
delta is defined as having a magnitude of infinity at an infinitesimally
small point $(x_0,y_0)$. The dirac spreads the scalar force at a time
$t$ in the function $f(t)$ to the appropriate point on the plate. The
condition for the dirac is given below: $$\begin{aligned}
\delta(x-x_0, y-y_0) =
  \begin{cases}
    \infty,       & \quad dS = 0\\
    0,  & \quad dS \neq 0\\
  \end{cases}, \quad \int_D{\delta(x-x_0, y-y_0)} = 1\end{aligned}$$

In discrete space, there is no longer an infinitesimally small point to
be considered. There is now a minimum grid spacing in the form of $h$.
In order for the dirac to have the same behaviour, its magnitude needs
to be changed to $\frac{1}{h^2}$. This is shown in Figure
\[figs:dirac\], where if $h_x = h_y = h$ then $h_xh_y = h^2$. In the FD
scheme this is implemented with a spreading function J. For the time
being in this case J is a flooring function \[eq:floorJ\]. Truncating
the force co-ordinates to the nearest grid point on the plate.

![The dirac delta in a) continuous space and b) discrete space.<span
data-label="figs:dirac"></span>](../Chapter_2/_Figs/diracdelta){width="14cm"}

$$\label{eq:floorJ}
  J_0=
  \begin{blockarray}{*{1}{c} l}
    \begin{block}{*{1}{>{$\footnotesize}c<{$}} l}
          \end{block}
    \begin{block}{[*{1}{c}]>{$\footnotesize}l<{$}}
      0  \bigstrut[t]&  \\
      \vdots &  \\
            1 & ($x_0,y_0$) \\
           \vdots &  \\
                      0 &  \\
    \end{block}
  \end{blockarray}$$

J can also be an interpolation function of some variety but this will be
discussed further in section \[chapter3:sec\_interp\]. For clarity the
$1/h^2$ is kept outside of J however, it can be included, along with
other coefficients when coded (see Appendix \[cd:plate\_string\_FD\]
lines 201-251).

$$\delta_{tt}u = -\kappa^{2}\delta_{\Delta\Delta}u^n- 2\sigma_{0}\delta_{t\cdot}u + \frac{1}{\rho H h^2}Jf(n)\nonumber \\$$

The force is now some discrete function of sample index $n$. A rough
modelling of a mallet would be a raised cosine over a period of time
(Figure \[figs:rc\]). The force does not need to simply be a simple
function, it can be any stream of numbers. Exchanging the raised cosine
for a `.wav` file will produce a plate reverb effect. Choosing two
readout points on the plate will produce a stereo plate reverb which can
then be mixed with the original audio. The plate model as instrument can
then be repurposed as a digital audio effect

![Force function f(n) as a raised cosine<span
data-label="figs:rc"></span>](../Chapter_2/_Figs/rc){width="10cm"}

### Modal Analysis

A second sanity check to ensure that the plate scheme is adhering to the
behaviour of the physical model is through modal analysis. If the modes
of vibration of FD model reflect the predicted modes we can safely
assume there are no inherent errors in the scheme. Under simply
supported conditions, prediction of the modes of vibration is reasonably
simplistic.

Firstly an ansatz of the plate is predicted as
$u(x,y,t) = e^{ik_{x}x}e^{ik_{y}y}e^{i\omega t}$ which is is separable
to $u(x,y,t)= T(t)X(x)Y(y)$. Since the boundary points can only ever be
equal to zero and the movement of the plate is to be oscillatory, it can
be assumed that behaviour across $x$ and $y$ can be represented as sine
waves.

$$\label{eq:ssmodes}
f_{p,q}  = \frac{\pi\kappa^{2}}{2}\left(\frac{p^{2}}{L_x} + \frac{q^{2}}{L_y}\right)$$

![Simply Supported Modal Frequencies at a) no oversampling and b) 8
times oversampling<span
data-label="figs:ss_modesfreqs"></span>](../Chapter_2/_Figs/SSOSR1-8){width="9cm"}

In \[eq:ssmodes\], $f$ is a frequency in $Hz$ and $p,q$ is the mode
number (see Appendix \[dvn:modes\]). To check that the scheme adheres to
these mode numbers, the frequency output can be compared to the
pre-calculated modes. Figure \[figs:ss\_modesfreqs\] shows the frequency
response of the scheme under two sampling rates. What we find is a close
match in lower frequencies. As we go up the frequency spectrum, there is
a drift between the solution and scheme output as a result of numerical
dispersion. This is rectified partially with increasing the sample rate,
Figure \[figs:ss\_modesfreqs\] (b), though computational expense is
increased greatly. Adjusting the grid spacing to take into account the
oversampling will undo any benefits from oversampling in the first
place. If using the plate model for musical purposes, the benefits of
oversampling will be minimal, especially if trying to keep within real
time computation. It is worth noting that the modes present will depend
on both the read-in and readout position of the plate; should either lie
on an anti-nodal point of frequency it will not be present in the
output.

![Clamped Modal Frequencies<span
data-label="figs:clamp_modes"></span>](../Chapter_2/_Figs/clamped-square-plate){width="8cm"}

![Thin Plate FD Modal Frequency Response and Approximated Modal
Frequencies under clamped conditions<span
data-label="figs:cl_modes"></span>](../Chapter_2/_Figs/ClampedModes){width="9cm"}

Under clamped conditions it becomes a little more difficult to analyse
the modal behaviour. Using [@LeissaA.W.1969Vop] as a reference we can
compare the predicted modes. The modes for a square plate from
[@LeissaA.W.1969Vop Table 4.24, Page 61] have been reproduced in figure
\[figs:clamp\_modes\]. Figure \[figs:cl\_modes\] shows the output of the
thin plate FTDT model with equal length and the predicted modal
frequencies from \[figs:clamp\_modes\].

When driven with a sine wave matching a modal frequency, the plate
exhibits the checker board pattern that we would expect from clamped and
fixed boundary conditions. Figure \[figs:33mode\] shows the 3-3 mode of
the plate. Further mode patterns are illustrated in Figure
\[figs:modePatterns\] replicated from [@BilbaoNSS Figure 11.3, page
309]. Which all in all, further demonstrates that at the very least, the
FTDT plate model is within the bounds of expected behaviour.

![$f_{3,3}$ Modal frequency<span
data-label="figs:33mode"></span>](../Chapter_2/_Figs/Mode3-3){width="3cm"}

![Modal frequency patterns under simply supported boundary conditions
[@BilbaoNSS Figure 11.3]<span
data-label="figs:modePatterns"></span>](../Chapter_2/_Figs/modePatterns){width="9cm"}

Coupling Plate and Strings {#chapter3}
==========================

This chapter will concern the problem of coupling the FDTD plate scheme
discussed in Chapter \[chapter2\]. Also, it will look into the concepts
in used order to achieve this goal; first deriving coupling conditions
through energy analysis and then implementing them within the FDTD
scheme. There are a number of constants/coefficients that are common
between the plate and the string. For this chapter, and the remainder of
this paper, all terms relating to just the plate will have a subscript
$_p$ and all those relating to the string, a subscript $s$. For
instance, the plate density $\rho_p$ and string denisty $\rho_s$.

For coupling a plate and a string an assumption is made that the end of
a string is rigidly connected to some point on a plate. Any two systems
can be coupled together. This chapter should hopefully provide the
groundwork for how to approach the coupling of two FDTD systems. For
deriving the force that a single string exerts on a plate, we begin with
equation \[eq:PFEOM\]. Little alteration needs to be made for coupling
multiple strings, though this will be discussed in greater depth in
\[chapter3:sec\_mult\_s\].

$$\begin{aligned}
\ddot{u} &=& -\kappa^{2}\Delta\Delta u + \frac{1}{\rho_{p}H}\delta(x - x_{e}, y - y_{e})f \nonumber \\
\label{eq:PFEOM} $$

The model of applying force in \[eq:PFEOM\] is, again, with a dirac
delta. It is identical to the application of force to the plate
discussed in section \[chapter2:sec\_force\]. The difference here is
that the means for the string driving the plate must be found. It is no
longer a function of sample index f(n), but a scalar force exerted on
the plate by the string.

1D Wave and Plate
-----------------

The first string model that will be coupled to the plate is the 1D wave
model

$$\begin{aligned}
\rho_s A\ddot{w} &= Tw^{\prime\prime} \label{eq:1DW1}  \\ 
\ddot{w} &= c^2w^{\prime\prime}, \quad c = \sqrt{\frac{T}{\rho_s A}} \label{eq:1DW2}\end{aligned}$$

Here, $w$ repsents tanversal displacement of the string, $\rho_s$ is the
density of the string, $A$ the cross sectional area and $T$ is the
tension. The constant $c$ is the wave speed in the string. Like the
plate it is advisable to return to \[eq:1DW2\] as greater complexity is
added. A translation to a FDTD model of \[eq:1DW2\] is:-

$$\begin{aligned}
\label{eq:1DWaveFDTD}
\delta_{tt}w = c^2\delta_{\eta\eta}w.\end{aligned}$$

In this instance, $\eta$ refers to the spatial dimension of the string.
The spatial dimensions of the string and the plate are completely
uncorrelated. This can make for some interesting combinations that would
be impossible to realise physically. A ten meter long string where the
ends are both coupled a mm apart on one meter squared plate would be
very easy to mock up in a FD scheme. To couple the plate and string
together we turn to energy analysis.

### Coupling conditions

Rather than ‘boundary conditions’ there are ‘coupling conditions’, in
this case a rigid coupling between the plate and the string. Assuming
all energy is contained in both the string and the plate, equation
\[eq:coup\_eng\], the change in energy being equal to the sum of the
energy of the string and the plate leaves us with a boundary term for
the string and a force term for the plate. Equating the units (speed of
force, speed of boundary) to each other gives the first interpretation
of how to couple the two systems.

$$\label{eq:coup_eng}
\frac{dE}{dt} = \mathcal{B}_p + \mathcal{B}_s + \int_{\mathcal{D}_s}{Jf\dot{u}}$$

In this instance we know that the plate boundaries $\mathcal{B}_p$ equal
0 and we are only connecting one end of the string to the plate at
$w_0$. Given this is a lossy system the change in energy is zero,
therefore all of this reduces to a term of the force and a boundary
condition for the string. For the 1D Wave Equation, the boundary term
is:

$$\label{eq:1DBCs}
 \mathcal{B}_s = T[\dot{w}w^\prime]_{0}^{N_{s}} = \dot{w}_{N_s}w^{\prime}_{N_s} - \dot{w}_{0}w^{\prime_{0}}$$

If we aim to only connect one end of the string to the plate then in
equation \[eq:1DBCs\], $N_s$ is the unconnected end of the string (.i.e
$\dot{w}_{N_s}w^{\prime}_{N_s} = 0$). As such, it is free to have
typical simply supported or clamped boundary conditions applied to it.
The benefit of only connecting one end is that there is control over the
boundary conditions and as such, some aesthetic possibilities are left
open. Taking all this into account we can express the remaining terms as
in equation \[eq:PSEB\].

$$\label{eq:PSEB}
T\dot{w}_{N_s}w^{\prime}_{N_s} - T\dot{w}_{0}w^{\prime}_{0} + \int_{\mathcal{D}_s}{\delta f\dot{u}} = 0$$

Converting this to FDTD, \[eq:PSEB\] transforms to:-

$$\begin{aligned}
h_p^2 \sum_{\mathcal{D}_p}{\frac{1}{h_p^2}Jf\delta_{t\cdot}u} &= T\delta_{t\cdot}w_{0}\delta_{\eta-}w_{0}\nonumber\\
fJ^T[\delta_{t\cdot}u] &= \delta_{t\cdot}w_{0}\delta_{\eta-}w_{0}
  \end{aligned}$$

The sum of the spreading vector $J$ leaves just the scalar force $f$ and
coupling point of the plate $J^T[\delta_{t\cdot}u]$. Equating the units
of the terms that are remaining we can impose conditions to satisfy
coupling the two systems. both $f$ and $T\delta_{\eta-}w_{0}$ are in
Newtons and the two velocity terms, $\delta_{t\cdot}w_0$ and
$J^T[\delta_{t\cdot}u]$ are of course in $ms^{-1}$. One interpretation
of these coupling conditions is given in \[eq:1Dcoup1\] along with a
FDTD equivalent \[eq:1Dcoup2\].

$$\begin{aligned}
\dot{w_{0}} = J^{T}[\dot{u}], \quad f = Tw_{0}^\prime \label{eq:1Dcoup1}\\
 \delta_{t\cdot}w_{0} = J^{T}[\delta_{t\cdot}u], \quad f = T\delta_{\eta-}w_0 \label{eq:1Dcoup2}
 \end{aligned}$$

The conditions in \[eq:1Dcoup2\] provide the starting point from which
we can couple the two systems. The schemes for both plate and the string
need to be expressed in terms of a centre time difference at the
relevant grid point. These can then be equated which will provide a term
for the force. The centre time difference schemes for the plate and
string are shown in equations \[eq:centFDp\] and \[eq:centFDs\].

$$\begin{aligned}
J^{T}[\delta_{t\cdot}u] &=& J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u + \frac{k}{2\rho_{p}Hh^{2}_{p}}Jf] \label{eq:centFDp} \\
\delta_{t\cdot}w_{0} &=& \frac{kc^{2}}{2h_{s}} (\delta_{x+}w_{0}- \frac{1}{T}f)  + \delta_{t-}w_{0} \label{eq:centFDs}\end{aligned}$$

Using \[eq:1Dcoup2\] we can equate both these terms and express just in
terms of the force.

$$\begin{aligned}
f  &=& \mathcal{M}(\frac{c^{2}}{h_{s}}\delta_{x+}w_{0} + \frac{2}{k}\delta_{t-}w_{0} - J^{T}[-\kappa^{2}\delta_{\Delta\Delta}u + \frac{2}{k}\delta_{t-}u]) \nonumber \\
 \label{eq:1DPf}\end{aligned}$$

Where $\mathcal{M}$ can be though of as the mass ratio between the plate
and the string and is defined as:
$$\mathcal{M} = \frac{\rho_{s}Ah_{s} \rho_{p}Hh^{2}_{p}}{\rho_{s}Ah_{s} + \rho_{p}Hh^{2}_{p}}$$

This definition of the force in \[eq:1DPf\] can then be fed back into
the string and the plate. From this, an update of the coupled scheme can
be derived just in terms of previously calculated points.

$$\begin{aligned}
w^{n+1}_0 &=& \lambda^{2}D_{x+}w_{0}  + 2w - w^{n-1} - \frac{k^{2}}{\rho A h_{s}}f
\label{eq:1DwF} $$

$$\begin{aligned}
J^{T}[u^{n+1}] &=& J^{T}[-\mu^{2}D_{\Delta\Delta}u  + 2u - u^{n-1}] + \frac{k^2}{\rho_{p} H h_p}f
\label{eq:PwF} $$

At this point, given the differences between the terms for force, plate
update and the string update, it is easy to get lost in the coefficient
chaos, especially when deriving this for the first time. What
\[eq:1DwF\] and \[eq:PwF\] have in common is the $k^2$ term in front of
the force term $f$. Taking the $k^2$ and integrating it into the force
term means that we can express the force purely in terms of the
coefficient matrix and a slightly modified version of the matrix $C$.

$$\begin{aligned}
F  &=& \mathcal{M}(\lambda^{2}D_{x+}w_{0} + 2w_{0} - 2w^{n-1}_{0} - J^{T}[-\mu^{2} D_{\Delta^{2}}u  + 2u - 2u^{n-1} ]) \nonumber \\
  F  &=& \mathcal{M}(B_{s(0,:)}w - J^T[B_{p}u] + 2C_{s(0,:)}w^{n-1} - J^T[2Cu^{n-1}])
\label{eq:recForce}\end{aligned}$$

The (0,:) subscript in \[eq:recForce\] is an elegant way of
cherry-picking the first row of the string matrices. We will see in a
later section, that combining the two schemes together will simplify
this problem. To a degree, \[eq:recForce\] remains true for the
inclusion of loss and stiffness of the string, both of which will be
discussed in future sections. Explicitly the update for the string and
the plate can be expressed as follows.

$$\begin{aligned}
w^{n+1}_0 &=& \lambda^{2}D_{x+}w_{0}  + 2w - w^{n-1} - \frac{1}{\rho_{s} A h_{s}}F\nonumber \\
J^{T}[u^{n+1}] &=& J^{T}[-\mu^{2}D_{\Delta\Delta}u  + 2u - u^{n-1}] + \frac{1}{\rho_{p} H h_p}F \nonumber \\
\label{eq:couplesimple}\end{aligned}$$

### 1D Wave /w Loss

![The string connected to the plate. The red line illustrates the point
of connection between the two systems<span
data-label="figs:String Plate"></span>](../Chapter_3/_Figs/StringPlate){width="12cm"}

The 1D wave equation is limited to just a generic loss term. This is
more useful than just a lossless system and is a good place to start
generating some output with discernible notes. It can also provide an
opportunity to begin planning on how to start expanding a model without
rewriting the underlying infrastructure.
$$(1 + k\sigma_0)w^{n+1} = (\lambda^{2} D_{\eta\eta} + I)w^n - (1 - k\sigma_0)w^{n-1}
 \label{1DwLoss}$$

The inclusion of loss in \[1DwLoss\] has some implication on the
derivation of the force term. The following should hopefully convey that
these changes are fairly trivial with alteration to existing code. The
derivation for \[1DwLoss\] can be found in section \[eq:1DFDloss\] which
is worth consulting to understand how to manipulate schemes that include
loss.

$$\begin{aligned}
 \left(\frac{1}{\rho_{s}Ah_{s}(1+k\sigma_{0s})} + \frac{1}{\rho_{p}Hh^{2}_{p}(1+k\sigma_{0p})}\right)f  &=& \frac{\frac{c^{2}}{h_{s}}\delta_{x+}w_{0} + \frac{2}{k}\delta_{t-}w_{0}}{(1+k\sigma_{0s})} - J^{T}\left[\frac{-\kappa^{2}\delta_{\Delta\Delta}u + \frac{2}{k}\delta_{t-}u}{(1+k\sigma_{0p})}\right] \nonumber \\
\label{eq:1DforceLoss}\end{aligned}$$

The generic loss is a centre time difference so can be brought over to
the LHS and factorised out (See Appendix derivation
\[dvn:1DforceLoss\]). The structure of the operation is mostly
unchanged. The mass ratio now contains the extra loss terms added by the
plate and the string

$$\mathcal{M} = \left(\frac{1}{\rho_{s}Ah_{s}(1+k\sigma_{0s})} + \frac{1}{\rho_{p}Hh^{2}_{p}(1+k\sigma_{0p})}\right)^{-1}$$

Following the same approach as the lossless case, the coupling force is
also the same as before for the current time step matrix $B$ and a
little alteration must be made for the $C$ matrix.

$$\begin{aligned}
F  &=&\mathcal{M}\left( \frac{ \lambda D_{x+}w_{0}  + 2w^{n}_{0} - 2w^{n-1}_{0} }{(1+k\sigma_{0s})} - J^{T}[\frac{-\kappa^{2}\delta_{\Delta\Delta}u + \frac{2}{k}\delta_{t-}u}{(1+k\sigma_{0p})}]\right) \nonumber \\
 F  &=& \mathcal{M} (B_{s(0,:)}w - J^T[B_{p}u] + (C_s-I_s)_{(0,:)}w^{n-1} - J^T[(C_p - I_p)u^{n-1}])\nonumber \\
\label{eq:recForce}\end{aligned}$$

It is worth noting at this point that all loss modelled in this paper
will be implemented explicitly. That is to say, the loss matrix for both
the string and the plate, $A$, is strictly diagonal and does not include
any spatial differences. As a result the inverse of $A$ can be
precomputed and stacked on top of the $B$ and $C$ matrices to minimise
computational costs.

### 1D Wave /w Force {#chapter3:sec_stringForce}

The exact same approach to force discussed in section
\[chapter2:sec\_force\] can be applied to the string.
$$\delta_{tt}w = c^2 \delta_{\eta\eta}w - 2\sigma_0 \delta{t\cdot}w + \frac{1}{\rho_s A h}\mathcal{J}_{f} f(n)$$

The operator $\mathcal{J}_{f}$ just locates the force at right grid
point for the vector $w$. It can be either a raised cosine for strikes
or half raised cosine for plucks. The important point to remember when
applying any force is the $\frac{k^2}{\rho_s A h}$ coefficient that is
in front of any force term for FDTD string model. Also, there is nothing
to stop force from being applied to both plate and string for percussive
accompaniment. Imagination is the limit on this front.

![Force Vector a) overwritten b) summation: What should be an
imperceptible double strike (b) can result in an unpleasant pluck
(a).<span
data-label="figs:writeForce"></span>](../Chapter_3/_Figs/writeForce){width="12cm"}

The string force can also be interpreted as a series of vectors for each
string, containing raised cosines at the required excitation points. One
possible fault in this approach is the situation where two excitations
on the same string overlap. It is important that the excitation is not
overwritten as (figure \[figs:writeForce\] a.), as this will create a
sudden pluck. This will result in a very prominent note playing. An
unusual case but, worth being aware of especially if excitation is
generated by some data source such as a MIDI file. Figure
\[figs:writeForce\] illustrates the product of both approaches.

### Joining vectors

What can be beneficial for computation is to amalgamate the update for
the string and the plate into one neat operation. With a force acting on
the string and the string driving the plate it is possible to confine
all this into one operation
$A\vec{uw} = B\vec{uw} + C\vec{uw} + J_{f}\mathcal{F}_{n}$. where $J_f$
spreads a force across the correct point on the string. $\mathcal{F}$ is
a pre-computed vector which is accessed at sample index ‘n’. For the
case of multiple string $I$ is a matrix $N \times N_s$ and $\mathcal{F}$
is a matrix $N_f \times N_s$. $N$ is the total number of grid points
across the plate and all strings, $N_s$ is the number of strings and
$N_f$ is the number of samples to be computed. The case of multiple
strings will be discussed in section \[chapter3:sec\_mult\_s\].

The spreading operator, $J$, can be thought of as a concatenation of two
spreading vectors. One that spread a force onto the plate $J_p$ and
another that connects to the end of the string $J_s$. The equation for
the coupling force can be restated as:-

![The Matrix Pattern of a joined plate and string system. Coupled points
have been circled and the sections demarcated to reflect which system
they refer.<span
data-label="figs:CouplingSpy"></span>](../Chapter_3/_Figs/CouplingSpy){width="8cm"}

$$\begin{aligned}
 F  &=& \mathcal{M}(\underbrace{(J_{s}^{T}[B_{s}] - J_{p}^T[B_{p}])}_{F^{n}}\vec{uw}^n + \underbrace{(J_{p}^T[(C_p - I_p)]+ J_{s}^{T}(C_s-I_s))}_{F^{n-1}}\vec{uw}^{n-1}])\nonumber \\
 F^n &=& \mathcal{M}(J_{s}^{T}[B_{s}] - J_{p}^T[B_{p}])\nonumber \\
 F^{n-1} &=& \mathcal{M}(J_{p}^T[C_p - I_p]+ J_{s}^{T}[C_s-I_s])
\label{eq:recForceJoin}\end{aligned}$$

If the plate and the string are joined in a single vector $\vec{uw}$,
\[eq:recForceJoin\] can be integrated into the coefficient matrices for
the plate and string. In this case, $F^n$ and $F^{n-1}$ are the force
coefficients of the B and C matrices for coupling the plate and string.
Remembering that the concatenated vectors are now operated on by
concatenated matrices. If this is difficult to follow, see Appendix
\[cd:plate\_string\_FD\] lines 250 - 251. For spreading the force back
out, the coefficients in front of the force term in \[eq:couplesimple\]
can be integrated into the new spreading vector $J$ (see Appenix B
\[cd:plate\_string\_FD\] lines 216 - 238).

$$\begin{aligned}
\vec{uw}^{n+1} = \underbrace{(B_1 + JF^n)}_\text{\large{B}}\vec{uw}^{n}- \underbrace{(C_1 + JF^{n-1})}_\text{\large{C}}\vec{uw}^{n-1}+ J_{f}\mathcal{F}_{n}
\label{eq:ForceWMatrix}\end{aligned}$$

The Matrices $B_1$ and $C_1$ refer to the combined coefficient matrices
without the coupling force included. What \[eq:ForceWMatrix\] implies is
that the full Matrices $B$ and $C$ can be precomputed along with the
loss that comes from the $A$ matrix. It also allows for easy expansion
to multiple strings as well as beginning to include a more sophisticated
interpolation for connecting to the plate.

Once the force is integrated into the coefficient matrix the coupling
points should be visible. This is illustrated in figure
\[figs:CouplingSpy\], which shows the region of the matrix that relates
to the plate and the string and sparse areas that lie to the side
indicate that a point of one system is coupled to another.

### Interpolation {#chapter3:sec_interp}

![The cross represents the coupling point lieing between grid points.
The normalised distance (or percentage) that the coupling point lies
towards the next grid point on the x and y axis is represented by
$\alpha_x$ and $\alpha_y$ respectively.<span
data-label="figs:linInterp"></span>](../Chapter_3/_Figs/linInterp){width="10cm"}

It is worth mentioning a different approach to the $J$ operator at this
stage. A more flexible way of coupling the strings is by a simple linear
interpolation model. When selecting the desired coordinates of the
string coupling point, it may not lie directly on a grid point. Since
nothing technically exists between grid points, a method for finding a
likely value for four grid points must be found. The force must be
spread amongst 4 grid points proportionally (see Figure
\[figs:linInterp\]). $J$ is no longer just a 1 at a chosen co-ordinate
but reflects the normalised distance the coupling point lies between
each grid point, \[eq:linJ\].

$$\label{eq:linJ}
  J=
  \begin{blockarray}{*{1}{c} l}
    \begin{block}{*{1}{>{$\footnotesize}c<{$}} l}
          \end{block}
    \begin{block}{[*{1}{c}]>{$\footnotesize}l<{$}}
      \vdots \bigstrut[t] &  \\
      (1-\alpha_x)(1-\alpha_y) & ($x_0,y_0$) \\
      (1-\alpha_x)\alpha_y & ($x_0,y_1$) \\
      \vdots &  \\
      \alpha_x(1-\alpha_y) & ($x_1,y_0$) \\
      \alpha_{x}\alpha_y & ($x_1,y_1$) \\
    \vdots &  \\
    \end{block}
  \end{blockarray}$$

Multiple Strings {#chapter3:sec_mult_s}
----------------

Just one string coupled to a plate has limited musical application.
Multiple strings allows for a more diverse set of tonal possibilities.
Adding extra strings is also a quick win in that it gives a sympathetic
vibration effect to the output which can be put to further creative use
(see \[chapter4:sec\_stringVerb\]).

When considering multiple strings the coupling force $F$ is no longer a
scalar but a vector. The spreading operator J becomes a matrix and the
string coefficients also become a matrix. Sticking with the flooring
operator, (i.e. no interpolation), some issues can arise issue,
particularly if it is preferable to avoid involved dealings with
interpolation. Any two strings attached to the same coupling point will
result in an unstable model. Extra error checking to ensure stability
would need to be implemented. This would mean limiting the number of
strings, ensuring all coupling points are unique or extending the number
of grid points to accommodate all the strings (This is the approach
taken in `thin_plate_stiff_string_loss.m` of the digital archive.). At
which stage, it is beneficial to just to implement interpolation
(discussed in section \[chapter3:sec\_interp\]).

$$\begin{aligned}
 \underbrace{(\frac{1}{\rho_{s}Ah_{s}(1+k\sigma_{0s})}I + \frac{1}{\rho_{p}Hh^{2}_{p}(1+k\sigma_{0p})}J^{T}_{p}J_{p})}_{\mathcal{M}^{-1}}f  &=& \frac{c^{2}}{h_{s}}\delta_{x+}w_{0} + \frac{2}{k}\delta_{t-}w_{0} - ...\nonumber \\
 &&J^{T}[-\kappa^{2}\delta_{\Delta\Delta}u + \frac{2}{k}\delta_{t-}u] \nonumber \\
\label{eq:multi1DforceLoss}
 \end{aligned}$$

There is not a great deal of change between \[eq:1DforceLoss\] and
\[eq:multi1DforceLoss\]. The main difference here is that the
multiplication of $J^{T}_{p}J_{p}$ is not a scalar. Should all
connecting strings share no grid points in common then $J^{T}_{p}J_{p}$
is an identity matrix for $f$. Otherwise, $J^{T}_{p}J_{p}$ is a more
complex matrix that will need to be inverted. The trem
$\frac{1}{\rho_{s}Ah_{s}(1+k\sigma_{0s})}$ is also a diagonal matrix
made up of the corresponding parameters for each string. Should all
strings share the same $\frac{1}{\rho_{s}Ah_{s}(1+k\sigma_{0s})} $ it
can be considered a scalar. The mass ratio term $\mathcal{M}$ is now an
inverse matrix from $\mathcal{M}^{-1}$ . There are a number of means for
determining this. MATLAB’s ‘backlash operator’, or inverse matrix method
or a Jacobi method approach. The luxury in this case is that it is
static, so it does not need to be recomputed every sample. If this were
not the case, it would be worth considering the Jacobi method approach.
Other than the construction of $\mathcal{M}$, the Force is identical to
\[eq:recForceJoin\]. The greatest difficulty at this point is indexing,
which is more a question of tedium rather than any great mental
challenge. Careful collation of coupling indices should allow the
construction of the spreading operators, and $\mathcal{M}$, to be a
fairly trivial procedure.

The Stiff String
----------------

Developing the scheme further, the 1D wave can be substituted for the a
stiff string model

\[eq:stiffstring\_loss\]
$$\rho_s A\ddot{w} = Tw^{\prime\prime} - E_{s}Iw^{\prime\prime\prime\prime}\\
\ddot{w} = c^{2}w^{\prime\prime} - \kappa_{s}^{2}w^{\prime\prime\prime\prime}
\label{eq:stiffFD}$$

The stiffness of the string brings with it a forth order term along with
Young’s Modulus for the string $E_s$ and the string moment of inertia
$I_s$.

$$w^{n+1} = c^{2}k^{2} \delta_{\eta\eta}w - \kappa^{2}k^{2} \delta_{\eta\eta\eta\eta} - 2\sigma_{0}\delta_{t\cdot}u + 2w - w^{n-1} + \frac{hk^2}{\rho A}Jf(n)
\label{eq:stiffFDu}$$

The update for the stiff string FD scheme is \[eq:stiffFDu\]. With the
stiffness comes a new bound on the grid spacing.

$$h_\eta \geq \sqrt{\frac{c^2k^2+\sqrt{c^4k^4+16\kappa_s^2k^2}}{2}}
\label{eq:stiffFDu}$$

Splitting equation \[eq:stiffFDu\] into coefficient matrices gives the
following:-

\[fd:stiffstring\_loss\_solv\] $$\begin{aligned}
\underbrace{(I+k\sigma_{0}I)}_{A_s}w^{n+1} =& \underbrace{(\lambda^2 D_{\eta\eta} + \frac{2k\sigma_{1}}{h^{2}}D_{\eta\eta} + 2I - \mu^{2} D_{\eta\eta\eta\eta})}_{B_s}w^{n} - ... \\
&\underbrace{(\frac{2k\sigma_{1}}{h^{2}}D_{\eta\eta} + (1-k\sigma_{0})I)}_{C_s}w^{n-1} + \frac{k^{2}}{\rho A h}Jf(t)\end{aligned}$$

Only this section would need to be altered in a FD script. The force
term for the concatenated system \[eq:ForceWMatrix\] goes unchanged,
making this a reasonably easy extension to implement.

### Coupling with Stiff String

With energy analysis a new set of coupling conditions are derived for
the stiff string. The change in energy of the stiff string,
$\frac{d\mathcal{E}}{dt}$, is defined as:-

$$\begin{aligned}
\frac{d\mathcal{E}}{dt} = T[\dot{w}w^{\prime}]_{0}^{N_\eta} - E_sI[\dot{w}w^{\prime\prime\prime}]_{0}^{N_\eta} + E_sI[\dot{w}^{\prime}w^{\prime\prime}]_{0}^{N_\eta}
 \end{aligned}$$

Comparing the units of force term from the plate energy, the following
FDTD coupling conditions are derived.

$$\begin{aligned}
\delta{t\cdot}w_0 &= J^T[\delta_{t\cdot}u]\label{eq:stiffCoup1}\\
\delta_{\eta\eta}w_0 &= 0 \label{eq:stiffCoup2}\\
f &= T\delta_{\eta-}w_0 - EI\delta_{\eta-}\delta_{\eta\eta}w_0\label{eq:stiffCoup3}\end{aligned}$$

In order to maintain the same $A = B + C$ structure of the finite
difference update, construction of the matrices for the string need to
be altered. In order for the scheme to comply to \[eq:stiffCoup1\] then
the $D_{\eta\eta}$ should have the first row as all zero.

$$\begin{aligned}
D_{\eta\eta} & = \frac{1}{h^{2}} \begin{bmatrix}
               0 & 0 & &   \\[0.3em]
               1 & -2 & 1 & \\[0.3em]
                  & \ddots & \ddots & \ddots  \\[0.3em]
                  & & b & -2
             \end{bmatrix}
             \label{eq:boundmat}\end{aligned}$$

The fourth order matrix should be constructed by squaring the second
order. This should help to keep the construction of coefficient matrices
automated.

$$D_{\eta\eta\eta\eta} = D_{\eta\eta}D_{\eta\eta}$$

The force term also changes but not by a great degree.

$$\begin{aligned}
\left(\frac{1}{\rho_{s}Ah_{s}(1+k\sigma_{0s})}\bar{I} + \frac{1}{\rho_{p}Hh_{p}^{2}(1+k\sigma_{0p})}J^{t}J\right)f = (\frac{\frac{c^2}{h}\delta_{\eta+}w_{0} - \frac{\kappa^{2}}{h_{s}^{2}}\delta_{\eta\eta}w_{1} + \frac{k}{2}w_{0}}{1+k\sigma_{0s}} -...\nonumber \\ J^{T}\left[\frac{-\kappa^{2}\delta_{\Delta\Delta}u + 2\sigma\delta_{t-}\delta_{\Delta}u}{1+k\sigma_{0p} }\right])\nonumber \\\end{aligned}$$

Should the string matrices be prepared correctly then
\[eq:recForceJoin\] and \[eq:ForceWMatrix\] should still hold true. This
means that again, very little intervention is required to actually
implement the stiff string (see Appendix \[cd:plate\_string\_FD\] lines
180-200).

### Frequency Dependant Loss {#chapter3:sec_freq_loss}

With the stiffness term in the string, frequency dependent loss term can
also be added. A new term appears in \[eq:stiffFD\],
$2\sigma_{1}\dot{w}^{\prime\prime}$ which requires a different approach
in the definition of $\sigma_0$.

$$w^{n+1} = c^{2}k^{2} \delta_{\eta\eta}w - \kappa^{2}k^{2} \delta_{\eta\eta\eta\eta} - 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{t-}\delta_{\eta\eta} + 2w - w^{n-1} + \frac{hk^2}{\rho A}Jf(n)
\label{eq:stiffFDu}$$

The calculation for $\sigma_0$ now changes. Taking the same approach as
[@BilbaoNSS page 178] where $\sigma$ is a function of frequency. We
first derive a characteristic equation for the string using an ansatz:

$$w = e^{st} e^{j \beta \eta}, \quad s = \sigma + j \omega\\$$

from which we derive a function of frequency

$$\sigma(\omega) = -\sigma_0-\sigma_1 \xi(\omega), \quad\xi(\omega) = \frac{-c^2 + \sqrt{c^4 + \kappa^2\omega^2}}{2\kappa^2}$$

Setting a T$_60$ for 2 frequency to the same as what is given in
equation (7.29) of [@BilbaoNSS page 178].

$$\sigma_0  = \frac{6\ln(10)}{\xi(\omega_2)-\xi(\omega_1)}\left(\frac{\xi(\omega_2)}{T_{60}(\omega_1)} - \frac{\xi(\omega_2)}{T_{60}(\omega2)}\right), \quad \sigma_1  = \frac{6\ln(10)}{\xi(\omega_2)-\xi(\omega_1)}\left(-\frac{1}{T_{60}(\omega_1)} + \frac{1}{T_{60}(\omega2)}\right)
\label{eq:freqDepLoss}$$

For the plate, the exact same apporach can be taken except, there is no
$c^2$ term in the plate. As such, the $\xi$ function reduces to

$$\xi(\omega) = \frac{\omega}{\kappa}$$

$$u^{n+1}\underbrace{(I+k\sigma_{0}I)}_\text{\Large{A}} = \underbrace{(-\mu_p^{2}D_{\Delta\Delta}+2k\sigma_{1p}D_{\Delta} + 2I)}_\text{\Large{B}}u^{n}-\underbrace{((1-k\sigma_0)I + 2k\sigma_{1p}D_{\Delta})}_\text{\Large{C}}u^{n-1}$$

If $\sigma_1 = 0$, the scheme reduces down to the terms of the generic
loss case. A conditional statement for loss should allow for easy
switching between the two methods of loss. Example code is available at
Appendix \[cd:plate\_string\_FD\] lines 100-130.

Alteration needs to take place to the stability condition of the plate
and string when loss is introduced. For plate stability conditions with
frequency dependant loss, consult [@ATorin_PhD]. For the plate and
string the new minimum grid spacing have been replicated below.

$$\begin{aligned}
h_p \geq \sqrt{4\sigma_{1p}k + \sqrt{(4\sigma_{1p}k)^2 + 12\kappa_{p}^{2}k^2}}\label{eq:freqStability1}\\ 
h_s \geq \sqrt{c^2k^2 + 4\sigma_{1s}k + \sqrt{(c^2k^2 + 4\sigma_{1s}k)^{2} + 16\kappa_{s}^{2}k^{2}}}\label{eq:freqStability2}\end{aligned}$$

Scheme Output
-------------

Selecting a point on the plate and assigning it to an output vector will
suffice to produce an output. Amplitude output such as this will be very
low passed. Amplitude output tends to exhibit a DC-like off-set. It can
cause difficulties in avoiding clipping on output.

It is advisable to use velocity as the output; this helps to suppress
the low frequency bias of amplitude output. Of course, this system is
not radiating in a room, coupled to air, so any output directly from the
plate will just be an approximation.

Should an excitation be close to a coupling point, there can be very
large spikes in output vector. This will very much depend on grid
spacing, magnitude of force and plate dimensions. It is worth keeping in
mind, perhaps built into error checking, to accommodate or anticipate
the proximity of the excitation to the output.

Consideration should also be given to the frequency content of the
output, which will depend on where the output ’pickup’ is placed. If the
read-out point lies directly at a nodal or anti-nodal point then
frequencies can be entirely absent or overly prominent respectively.
Treating the read-out point as a normalised co-ordinate, careful coding
can help to avoid this kind of problem, regardless of plate dimensions.

Simply selecting two opposing points on a plate will produce a stereo
output of fairly good quality. Experimenting with multiple points that
are spread across the stereo field with equal power panning can also
produce some interesting audio, see [@RoadsCurtis1995Tcmt 457-461].

Extensions to the FD Plate and Coupled String {#chapter4}
=============================================

Though many of the main issues concerning implementation have been
covered in the previous chapters, there are still some extensions and
methods of approach to these models that will improve the quality of the
output and code. This chapter will first deal with methods for
optimisation as well as any helpful tips for those unfamiliar with
programming. With regard to extending the schemes, some extensions have
been carried out and are contained in the digital archive accompanying
this paper. This includes generating notes by parsing MIDI data which
was briefly touched on in section \[chapter3:sec\_stringForce\]. MIDI
offers the easiest, ready-made solution but, for a greater degree of
nuance, a music notation based markup language may also be considered.
There are also some further ideas which have not been tackled, which are
suggested in section \[chapter4:sec\_suggestions\].

Optimisation
------------

For those taking on a coding project such as this for the first time it
can be beneficial to take into consideration a few helpful tips before
coding. It is best to attempt to compartmentalise each section of the
scheme; isolating each part in this way should allow for easy alteration
without having to tread over old ground. This is obvious for seasoned
coders but still worth remembering, even for a fairly short script. Of
course there are always times when optimising will require amalgamating
or discarding entire sections and there is no real solution but to
rebuild the infrastructure of a script.

Taking advantage of some file versioning system is highly recommended.
Not keeping track of changes to code can have disastrous consequences.
Simply creating new files for every edit will also lead to a cluttered
library of scripts fairly quickly, though is better than no versioning
at all. Creating bespoke functions will greatly reduce space and code
duplication; some 2D FDTD functions have been suggested in Appendix
\[appB:sec1\]. A little tweaking should allow for an arbitrary number of
dimensions as well as boundary conditions. The best means of optimising
code in MATLAB is to keep coefficient matrices in sparse matrix format.
Careful linear indexing and a little thought can render it a fairly
simple problem (see section \[chapter1:sec\_LinDex\]). This is the
approach used for all FDTD schemes in the digital archive. Appendix
\[cd:plate\_string\_FD\] presents a simple single string/plate model.
From experience, it is tricky to initially get to grips with this kind
of logic. Starting with small and theoretical examples should help
familiarisation.

Joining systems together will greatly improve computation speed. This
has been discussed in the previous chapter, though it is worth
reiterating that the real difficulty here is careful indexing. The key
goal in any optimisation is to scrutinise every aspect of the model and
avoid duplicating any computation.

Calibration {#chapter4:sec_calibration}
-----------

An easy method for producing timbres that are less synthetic involves
calibrating the models against some real recordings. This has been
carried the case of the Piano Models in the digital archive. Individual
strings of a piano were recorded while all others were damped with
fabric. The frequency spectrum of the piano strings was compared to the
output of isolated FDTD stiff string. The stiff string parameters, such
as Young’s Modulus, length, radius, loss were altered until the
frequency response of the FDTD model closely resembled that of the
recording. Changes to Young’s modulus accommodate in part for changes in
string behaviour from overwinding the solid core with another material
(see Figure \[figs:Overwound\]).

![Overwound Piano String<span
data-label="figs:Overwound"></span>](../Chapter_4/_Figs/Overwound){width="8cm"}

It is advisable to construct an environment that will allow for easy
comparison between a string model and sound file. With sensible file
naming, a single changed variable should be all that is required to
compare different strings. The data can either be stored in a separate
file or supplant the current string parameters. One warning when
undertaking this is that the stability condition for grid spacing should
be changed, as in \[eq:freqStability1\] and \[eq:freqStability2\].
Without this change it is likely that the model will be unstable unless
under very particular conditions.

With regard to the plate properties, [@1999Wh:w Chapter 5] has been a
valuable reference of poisson ratios and densities of wood as well as
[@bremaud2011anisotropy] which is a useful source of shear coefficients.

Reverberation {#chapter4:sec_stringVerb}
-------------

The sympathetic vibration of multiple strings on a plate can lead to
some interesting reverberation effects. The effect of a piano with a
sustain pedal held down during recording was a favourite of Frank Zappa
[@carr2016frank 79] and driving the plate/string model with audio is one
approach to mimicking the behaviour of systems such as the Ondes
Martenot’s Palme Diffuseur (Figure \[figs:Overwound\]).

![The Palme Diffuseur for the Ondes Martenot [@ondes]<span
data-label="figs:Overwound"></span>](../Chapter_4/_Figs/Ondes.jpg){width="8cm"}

The FDTD scheme suggested may not be constrained to behaving as musical
instrument but may also be used as a digital audio effect. In this
paper, just one plate and one set of strings has been coupled. Their is
no real limit to the number of coupled systems and a whole realm of
physically impossible instruments and effects is available.

Further additions {#chapter4:sec_suggestions}
-----------------

Aesthetically the output of the previously discussed schemes is still
noticeably synthetic. This is not to say that there are not musical
possibilities however, if the goal of the reader is to create sounds
approaching a ‘real’ instrument there are a number of ways forward. For
fretted instruments there are the interactions between the string and
the fretboard [@bilbao2015numerical]. Damping, mallet and other force
modelling is covered well in [@BilbaoNSS], as well as nonlinear models
of both the string and plate. A variation in the string model will
provide the greatest aesthetic change, so may be the first place to
consider investigating.

For the plate model, there is the addition of anisotropy. The difficulty
here is the decision of the shear modulus for a material, of which there
is no consensus. This becomes a more intense problem to solve and may
not be for the casual physical modeller. Coupling with air and irregular
geometries will require a step towards finite element or finite volume
approaches both on which are beyond the scope of this paper.

Reflections
-----------

The initial goal of this project was to create a full irregularly-shaped
air-coupled instrument model. Under the constraints of a three month
working period, this has proven to be a little too ambitious. With what
has been created instead, there is still potential for musical
applications. Throughout the project the main division has been between
analysis and implementation. With a limited time scope, to implement the
most diverse and interesting model requires that it be checked for
stability, and validity, analytically. To maximise interesting results,
a healthy balance has to be found between the two. Hopefully this paper
will have provided some use to those undertaking a similar project, and
helped to illuminate the more difficult concepts more clearly.

Full Derivation of Formulas {#appA}
===========================

Each section of the Appendix relates to the equivalent Chapter

The Thin Plate {#appA:sec2}
--------------

### Lossless Thin Plate

$$\begin{aligned}
\rho H \delta_{tt}u &= -D\delta_{\Delta\Delta}u \nonumber \\
\delta_{tt}u &= -\kappa^{2}\delta_{\Delta\Delta}u\nonumber \\
\frac{1}{k^2}(u^{n+1} - 2u^{n} + u^{n-1}) &= -\kappa^{2}\delta_{\Delta\Delta}u\nonumber \\
u^{n+1} &= -k^{2}\kappa^{2}\delta_{\Delta\Delta}u + 2u^{n} - u^{n-1} \nonumber \\
 &= -\mu^2 D_{\Delta\Delta}u + 2u^{n} - u^{n-1} \nonumber \\\end{aligned}$$

### Thin Plate w/ Generic Loss {#dvn:plateFDTDWloss}

$$\begin{aligned}
\rho H \delta_{tt}u &= -D\delta_{\Delta\Delta}u - 2\rho H \sigma_{0}\delta_{t\cdot}u \nonumber \\
\delta_{tt}u &= -\kappa^{2}\delta_{\Delta\Delta}u^n- 2\sigma_{0}\delta_{t\cdot}u\nonumber \\
\frac{1}{k^2}(u^{n+1} - 2u^{n} + u^{n-1}) &= -\kappa^{2}\delta_{\Delta\Delta}u^n-\frac{2\sigma_{0}}{2k}u^{n+1} + \frac{2\sigma_{0}}{2k}u^{n-1}\nonumber \\
u^{n+1} &= -k^{2}\kappa^{2}\delta_{\Delta\Delta}u^n  + k^{2}\sigma_{0}u^{n-1} + 2u^{n} - u^{n-1} \nonumber \\
u^{n+1}+ k^{2}\sigma_{0}u^{n+1} &= -k^{2}\kappa^{2}\delta_{\Delta\Delta}u^n + k^{2}\sigma_{0}u^{n-1} + 2u^{n} - u^{n-1} \nonumber \\
\underbrace{(1+ k^{2}\sigma_{0})}_{A}u^{n+1} &= \underbrace{(-\mu^2 D_{\Delta\Delta} + 2I)}_{B}u^{n} - \underbrace{(1 - k^{2}\sigma_{0})}_{C}u^{n-1} \nonumber \\\end{aligned}$$

### Thin Plate w/ Frequency Dependant Loss

$$\begin{aligned}
\rho H \ddot{u}  &=& -D\Delta\Delta u -2\rho H \sigma_{0}\dot{u} +2\rho H \sigma_{1}\Delta\dot{u}\nonumber \\
\ddot{u} &=& -\kappa^{2}\Delta\Delta u -2\sigma_{0}\dot{u} +2\sigma_{1}\Delta\dot{u}\nonumber \\
\delta_{tt}u &=& -\kappa^{2}\delta_{\Delta\Delta}u +- 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{\Delta}\delta_{t-}u \nonumber \\
\frac{1}{k^2}(u^{n+1} - 2u^{n} + u^{n-1}) &=& -\kappa^{2}\delta_{\Delta\Delta}u - 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{\Delta}\delta_{t-}u\nonumber \\
u^{n+1} &=& k^2(-\kappa^{2}\delta_{\Delta\Delta}u - 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{\Delta}\delta_{t-}u +  \frac{1}{\rho_{p}Hh^{2}_{p}}Jf) + 2u^{n} - u^{n-1} \nonumber \\
u^{n+1} &=& -\kappa^{2}k^2\delta_{\Delta\Delta}u - \frac{2\sigma_{0}k^2}{2k}(u^{n+1}-u^{n+1}) + ... \nonumber \\
&&2\sigma_{1}k^2\delta_{\Delta}\delta_{t-}u \\
u^{n+1} + \sigma_{0}ku^{n+1} &=& -\kappa^{2}k^2\delta_{\Delta\Delta}u + 2u^{n}  + 2\sigma_{1}k\delta_{\Delta}u + ... \nonumber \\
&& \sigma_{0}ku^{n-1}-  u^{n-1} -2\sigma_{1}k\delta_{\Delta}u^{n-1} \nonumber \\
\underbrace{(I+k\sigma_{0}I)}_\text{\large{A}}u^{n+1} &=& \underbrace{(-\mu^{2}D_{\Delta\Delta}+2\sigma_{1}D_{\Delta} + 2I)}_\text{\large{B}}u^{n}-...\nonumber \\
&&\underbrace{((1-k\sigma_0)I + 2\sigma_{1}D_{\Delta})}_\text{\large{C}}u^{n-1} \end{aligned}$$

### The Bihamonic Operator {#the-bihamonic-operator .unnumbered}

$$\begin{aligned}
\delta_{\Delta\Delta} & = & (\delta_{\Delta})^{2}\nonumber \\
(\delta_{\Delta})^{2} & = & (e_{x+1} - 2 + e_{x-1} e_{y+1} - 2 + e_{y-1})^{2} \nonumber \\ 
 & = & (20 - 8(e_{x+} + e_{x-} + e_{y+} + e_{y-})  + ... \nonumber \\
&& (e_{x+,y+} + e_{x+,y-} + e_{x-,y+} + e_{x-,y-}) + ... \nonumber \\
&& (e_{x+2} + e_{x-2} + e_{y+2} + e_{y-2})\end{aligned}$$

\[dvn:boundaries\] $$\begin{aligned}
\text{Simply Supported:} \nonumber \\
\delta_{xx}u_{0} &=& 0 \nonumber \\
\delta_{xx}u_{0} &=& \frac{1}{h^{2}}\left(u_{1} - 2u_{0} + u_{-1}\right) \nonumber \\
\frac{1}{h^{2}}\left(u_{1} - 2u_{0} + u_{-1}\right) &=& 0 \nonumber \\
u_{1} - 2u_{0} + u_{-1} &=& 0  \nonumber \\
u_{-1} &=& 2u_{0} - u_{1}, \quad (u_{0} = 0)    \nonumber \\
u_{-1} &=& - u_{1}\end{aligned}$$ $$\begin{aligned}
\text{Clamped:} \nonumber \\
\delta_{x\cdot}u_{0} &=& 0 \nonumber \\
\delta_{x\cdot}u_{0} &=& \frac{1}{2h}\left(u_{1} - u_{-1}\right) \nonumber \\
\frac{1}{2h}\left(u_{1} - u_{-1}\right)  &=& 0 \nonumber \\
u_{1} - u_{-1} &=& 0  \nonumber \\
u_{-1} &=& u_{1}\end{aligned}$$

### Energy: Continuous

Boundary conditions = $\mathcal{B}$ \[eq:PlateEnergy\] $$\begin{aligned}
\ddot{u} &=& -\kappa^{2}\Delta\Delta u \quad (\times \dot{u}, \int_D) \nonumber  \\
\int_{\mathcal{D}}{\dot{u}\ddot{u}} &=& -\kappa^{2}\int_{\mathcal{D}}{\dot{u}\Delta\Delta u} \nonumber  \\
\int_{\mathcal{D}}{\frac{\partial}{\partial t}[\frac{1}{2}\dot{u}^2]} =-\kappa^{2}\langle\dot{u},\Delta\Delta u\rangle_{\mathcal{D}} \nonumber  \\
\frac{\partial}{\partial t}\left[\frac{1}{2}\int_{\mathcal{D}}{\dot{u}^2}\right] &=&\mathcal{B} - \frac{\partial}{\partial t}\left[\frac{\kappa^{2}}{2}\langle\Delta u,\Delta u\rangle_{\mathcal{D}} \right] \nonumber  \\
\frac{\partial}{\partial t}\underbrace{\left[\frac{1}{2}\|\dot{u}\|^{2}_{\mathcal{D}} + \frac{\kappa^2}{2}\|\Delta u\|^{2}_{\mathcal{D}}\right]}_{E}  &=& \mathcal{B} 
$$

### Energy: Discrete {#dvn:energyDisc}

$$\begin{aligned}
\delta_{tt}u &=& -\kappa^{2}\delta_{\Delta\Delta} u \nonumber \\
h^2\sum_{\mathcal{D}}{\delta_{t\cdot}u\delta_{tt}u} &=& -\kappa^{2}h^2\sum_{\mathcal{D}}{\delta_{t\cdot}u\Delta\Delta u}\nonumber \\
h^2\sum_{\mathcal{D}}{\frac{1}{2}\delta_{t+}[\delta_{t-}u]^2} &=&-\kappa^{2}h^2\langle\delta_{t\cdot}u,\Delta\Delta u\rangle_{\mathcal{D}} \nonumber\\
\frac{h^2}{2}\delta_{t+}\|\delta_{t-}u\|^2 &=&-\kappa^{2}h^2\langle\delta_{t\cdot}u,\Delta\Delta u\rangle_{\mathcal{D}} \nonumber\\
\frac{h^2}{2}\delta_{t+}\|\delta_{t-}u\|^2&=&\mathcal{B} - \delta_{t+}\left[\frac{\kappa^{2}h^2}{2h^4}\langle D_{\Delta} u^n,D_{\Delta} u^{n-1}\rangle_{\mathcal{D}} \right]\nonumber\\
$$

### Energy: Discrete Solution {#energy-discrete-solution .unnumbered}

$$\delta_{t+}\left[\frac{h^2}{2k^2}(u^n-u^{n-1})^T(u^n-u^{n-1}) + \frac{\kappa^2}{2h^2}(D_{\Delta} u^n)^T(D_{\Delta} u^{n-1})\right] =0$$

### Modal Analysis: Simply Supported {#dvn:modes}

$$\begin{aligned}
u(x,y,t) &= e^{iwt} X(x)Y(y) \nonumber \\
X(0) &= X(L_{x})=0, \quad Y(0) = Y(L_{y}) = 0  \nonumber \\
X(x) &= \sin{(k_{x}x)}, \quad Y(y) = \sin{(k_{y}y)} \nonumber \\
k_{x} &= \frac{p\pi}{L_{x}}, \quad k_{y} = \frac{q\pi}{L_{y}}, \quad p,q = \mathbb{Z}^{+} \nonumber \\
\nonumber \\
\ddot{u} &= -\kappa^{2}\left[\pdv[4]{}{x} + 2\frac{\partial^{4}}{{\partial^{2}x}{\partial^{2}y}} +  \pdv[4]{}{y} \right]u \nonumber \\
-\omega^{2}u &= -\kappa^{2}\left(k_{x}^{4} + 2k_{x}^{2}k_{y}^{2} + k_{y}^{4} \right) u \nonumber \\
\omega^{2}u &= \kappa^{2}\left(k_{x}^{2} + k_{y}^{2} \right)^{2} u \nonumber \\
\omega &= \kappa \left(k_{x}^{2} + k_{y}^{2} \right) \nonumber \\
2\pi f_{p,q} &= \kappa \left(\frac{p\pi}{L_{x}}^{2} + \frac{q\pi}{L_{y}}^{2} \right) \nonumber \\
2\pi f_{p,q} &= \pi^{2}\kappa \left(\frac{p}{L_{x}}^{2} + \frac{q}{L_{y}}^{2} \right) \nonumber \\
f_{p,q}  &= \frac{\pi\kappa^{2}}{2}\left(\frac{p^{2}}{L_x} + \frac{q^{2}}{L_y}\right) \end{aligned}$$

### Lossless w/ Force {#plateFD}

$$\begin{aligned}
\ddot{u} &=& -\kappa^{2}\Delta\Delta u + \frac{1}{\rho_{p}H}\delta(x - x_{e}, y - y_{e})f \nonumber \\
\delta_{tt}u &=& -\kappa^{2}\delta_{\Delta\Delta}u + \frac{1}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
\frac{1}{k^{2}}(u^{n+1} - 2u + u^{n-1}) &=& -\kappa^{2}\delta_{\Delta\Delta}u + \frac{1}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
u^{n+1} &=& -k^{2}\kappa^{2}\delta_{\Delta\Delta}u + 2u - u^{n-1} + \frac{k^{2}}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
u^{n+1} &=& -\mu^{2}D_{\Delta\Delta}u + 2u - u^{n-1} + \frac{k^{2}}{\rho_{p}Hh^{2}_{p}}Jf\end{aligned}$$

### Frequency Dependant Loss w/ Force

$$\begin{aligned}
\rho H \ddot{u}  &=& -D\Delta\Delta u -2\rho H \sigma_{0}\dot{u} +2\rho H \sigma_{1}\Delta\dot{u} + \delta(x-x_{e},y-y_{e})f(t)\nonumber \\
\ddot{u} &=& -\kappa^{2}\Delta\Delta u -2\sigma_{0}\dot{u} +2\sigma_{1}\Delta\dot{u} + \frac{1}{\rho_{p}H}\delta(x - x_{e}, y - y_{e})f(t)\nonumber \\
\delta_{tt}u &=& -\kappa^{2}\delta_{\Delta\Delta}u +- 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{\Delta}\delta_{t-}u +  \frac{1}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
\frac{1}{k^2}(u^{n+1} - 2u^{n} + u^{n-1}) &=& -\kappa^{2}\delta_{\Delta\Delta}u - 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{\Delta}\delta_{t-}u +  \frac{1}{\rho_{p}Hh^{2}_{p}}Jf\nonumber \\
u^{n+1} &=& k^2(-\kappa^{2}\delta_{\Delta\Delta}u - 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{\Delta}\delta_{t-}u +  \frac{1}{\rho_{p}Hh^{2}_{p}}Jf) + 2u^{n} - u^{n-1} \nonumber \\
u^{n+1} &=& -\kappa^{2}k^2\delta_{\Delta\Delta}u - \frac{2\sigma_{0}k^2}{2k}(u^{n+1}-u^{n+1}) + ... \nonumber \\
&&2\sigma_{1}k^2\delta_{\Delta}\delta_{t-}u +  \frac{k^2}{\rho_{p}Hh^{2}_{p}}Jf + 2u^{n} - u^{n-1} \nonumber \\
u^{n+1} + \sigma_{0}ku^{n+1} &=& -\kappa^{2}k^2\delta_{\Delta\Delta}u + 2u^{n}  + 2\sigma_{1}k\delta_{\Delta}u + ... \nonumber \\
&& \sigma_{0}ku^{n-1}-  u^{n-1} -2\sigma_{1}k\delta_{\Delta}u^{n-1} + \frac{k^2}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
\underbrace{(I+k\sigma_{0}I)}_\text{\large{A}}u^{n+1} &=& \underbrace{(-\mu^{2}D_{\Delta\Delta}+2\sigma_{1}D_{\Delta} + 2I)}_\text{\large{B}}u^{n}-...\nonumber \\
&&\underbrace{((1-k\sigma_0)I + 2\sigma_{1}D_{\Delta})}_\text{\large{C}}u^{n-1} + \frac{k^2}{\rho_{p}Hh^{2}_{p}}Jf\end{aligned}$$

Coupling String and Plate {#appA:sec_coupling}
-------------------------

### FDTD 1D Wave Lossless

$$\begin{aligned}
\ddot{w} &=& c^2w^{\prime\prime} \nonumber \\
\delta_{tt}w &=& c^2 \delta_{\eta\eta}w \nonumber \\
\frac{1}{k^{2}}(w^{n+1} - 2w + w^{n-1}) &=& c^{2} \delta_{\eta\eta}w \nonumber \\
w^{n+1} &=& \lambda^{2} D_{\eta\eta}w  + 2w - w^{n-1}
\label{eq:1DFD}\end{aligned}$$

### FDTD 1D Wave Generic Loss

$$\begin{aligned}
\ddot{w} &=& c^2w^{\prime\prime} - 2\sigma_0 \dot{w} \nonumber \\
\delta_{tt}w &=& c^2 \delta_{\eta\eta}w -2\sigma_0 \delta{t\cdot}w \nonumber \\
\frac{1}{k^{2}}(w^{n+1} - 2w + w^{n-1}) &=& c^{2} \delta_{\eta\eta}w \nonumber \\
\underbrace{(1 + k\sigma_0)}_{\text{\large{A}}}w^{n+1} &=& \underbrace{(\lambda^{2} D_{\eta\eta} + I)}_{\text{\large{B}}}w^n - \underbrace{(1 - k\sigma_0)}_{\text{\large{C}}} w^{n-1}
\label{eq:1DFDloss}\end{aligned}$$

### Coupling Conditions: Energy Analysis

### 1D Wave Energy {#d-wave-energy .unnumbered}

$$\begin{aligned}
\int_{\mathcal{D}_s}{\rho_{s}A \dot{w}\ddot{w}} &=& \int_{\mathcal{D}_s}{T \dot{w}w^{\prime\prime}} \nonumber \\
\int_{\mathcal{D}_s}{\rho_{s}A \dot{w}\ddot{w}} &=&[\dot{w}w^{\prime}]_{0}^{N_s}- \int_{\mathcal{D}_s}{T \dot{w}^\prime w^{\prime}}\nonumber \\
\frac{dE_s}{dt} &=& [\dot{w}w^{\prime}]_{0}^{N_s} = \mathcal{B}_s\end{aligned}$$

### Plate and String Energy {#plate-and-string-energy .unnumbered}

$$\begin{aligned}
\frac{dE}{dt} &=& \mathcal{B}_p + \mathcal{B}_s + \int_{\mathcal{D}_s}{Jf\dot{u}} \nonumber \\
\int_{\mathcal{D}_s}{Jf\dot{u}} &=& fJ^T\dot{u}\end{aligned}$$

### Plate Centre Difference

$$\begin{aligned}
J^{T}\delta_{tt}u &=& J^{T}[-\kappa^{2}\delta_{\Delta\Delta}u + \frac{1}{\rho_{p}Hh^{2}_{p}}Jf] \nonumber \\
J^{T}[\frac{2}{k}(\delta_{t\cdot}u - \delta_{t-}u)] &=& J^{T}[-\kappa^{2}\delta_{\Delta\Delta}u + \frac{1}{\rho_{p}Hh^{2}_{p}}Jf] \nonumber \\
J^{T}[\delta_{t\cdot}u] &=& J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u + \frac{k}{2\rho_{p}Hh^{2}_{p}}Jf] \nonumber \\
\label{eq:platecentFD}\end{aligned}$$

### 1D Wave Centre Difference

$$\begin{aligned}
\delta_{tt}w &=& c^2 \delta_{xx}w \nonumber \\
\frac{2}{k}(\delta_{t\cdot}w_{0} - \delta_{t-}w_{0}) &=& c^2 \delta_{xx}w_{0} \nonumber \\
\delta_{t\cdot}w_{0} &=& \frac{kc^{2}}{2} \delta_{xx}w_{0}  + \delta_{t-}w_{0} \nonumber \\
&=& \frac{kc^{2}}{2h_{s}} (\delta_{x+}w_{0}- \frac{1}{T}f)  + \delta_{t-}w_{0}
\label{stringcentreFD}\end{aligned}$$

### 1D Wave/Plate Coupling Force

$$\begin{aligned}
 \frac{kc^{2}}{2h_{s}} (\delta_{x+}w_{0}- \frac{1}{T}f)  + \delta_{t-}w_{0} &=& J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u + \frac{k}{2\rho_{p}Hh^{2}_{p}}Jf] \nonumber \\
 (\frac{k}{2\rho_{p}Hh^{2}_{p}} + \frac{kc^{2}}{2Th_{s}})f  &=& \frac{kc^{2}}{2h_{s}}\delta_{x+}w_{0} + \delta_{t-}w_{0} - J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u] \nonumber \\
 \frac{k}{2}(\frac{1}{\rho_{s}Ah_{s}} + \frac{1}{\rho_{p}Hh^{2}_{p}})f  &=& \frac{kc^{2}}{2h_{s}}\delta_{x+}w_{0} + \delta_{t-}w_{0} - J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u] \nonumber \\
  (\frac{1}{\rho_{s}Ah_{s}} + \frac{1}{\rho_{p}Hh^{2}_{p}})f  &=& \frac{c^{2}}{h_{s}}\delta_{x+}w_{0} + \frac{2}{k}\delta_{t-}w_{0} - J^{T}[-\kappa^{2}\delta_{\Delta\Delta}u + \frac{2}{k}\delta_{t-}u] \nonumber \\
f &=& (\mathcal{M})(\frac{c^{2}}{h_{s}}\delta_{x+}w_{0} + \frac{2}{k}\delta_{t-}w_{0} - ...\nonumber \\ && J^{T}[-\kappa^{2} \delta_{\Delta^{2}}u  + \frac{2}{k}\delta_{t-}u ]))\nonumber \\
\label{dvn:forcedef}\end{aligned}$$

### Thin Plate w/ Coupling Force

$$\begin{aligned}
\ddot{u} &=& -\kappa^{2}\Delta\Delta u + \frac{1}{\rho_{p}H}\delta(x - x_{e}, y - y_{e})f \nonumber \\
\delta_{tt}u &=& -\kappa^{2}\delta_{\Delta\Delta}u + \frac{1}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
\frac{1}{k^{2}}(u^{n+1} - 2u + u^{n-1}) &=& -\kappa^{2}\delta_{\Delta\Delta}u + \frac{1}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
u^{n+1} &=& -k^{2}\kappa^{2}\delta_{\Delta\Delta}u + 2u - u^{n-1} + \frac{k^{2}}{\rho_{p}Hh^{2}_{p}}Jf \nonumber \\
 &=& -\mu^{2}D_{\Delta\Delta}u + 2u - u^{n-1} + \frac{k^{2}}{\rho_{p}Hh^{2}_{p}}Jf
\label{plateFD}\end{aligned}$$

### 1D Wave FD Scheme w/ Coupling Force

$$\begin{aligned}
w^{n+1}_0 &=& \frac{k^{2}c^{2}}{h_{s}} (\delta_{x+}w_0 - \frac{1}{T}f)  + 2w - w^{n-1} \nonumber \\
 &=& \lambda^{2}D_{x+}w_{0}  + 2w - w^{n-1} - \frac{k^{2}}{\rho A h_{s}}f\nonumber \\
 &=& \lambda^{2}D_{x+}w_{0}  + 2w - w^{n-1} - ... \nonumber \\ && \frac{\mathcal{M}}{\rho_{s}Ah_{s}} (\lambda^2\delta_{x+}w_{0} + 2w_{0} - 2w^{n-1}_{0} - ...\nonumber \\ && J^{T}[-\mu^{2} D_{\Delta^{2}}u  + 2u - 2u^{n-1} ])\nonumber \\
 &=& \lambda^{2}D_{x+}w_{0}  + 2w - w^{n-1} - ... \nonumber \\ && \mathcal{M}_{s}(\lambda^{2}D_{x+}w_{0} + 2w_{0} - 2w^{n-1}_{0} - ...\nonumber \\ && J^{T}[-\mu^{2} D_{\Delta^{2}}u  + 2u - 2u^{n-1} ]))\nonumber \\
\label{stringforce}\end{aligned}$$

### 1D Wave/Plate Coupling Force

$$\begin{aligned}
 \frac{kc^{2}}{2h_{s}} (\delta_{x+}w_{0}- \frac{1}{T}f)  + \delta_{t-}w_{0} &=& J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u + \frac{k}{2\rho_{p}Hh^{2}_{p}}Jf] \nonumber \\
 (\frac{k}{2\rho_{p}Hh^{2}_{p}} + \frac{kc^{2}}{2Th_{s}})f  &=& \frac{kc^{2}}{2h_{s}}\delta_{x+}w_{0} + \delta_{t-}w_{0} - J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u] \nonumber \\
 \frac{k}{2}(\frac{1}{\rho_{s}Ah_{s}} + \frac{1}{\rho_{p}Hh^{2}_{p}})f  &=& \frac{kc^{2}}{2h_{s}}\delta_{x+}w_{0} + \delta_{t-}w_{0} - J^{T}[-\frac{k\kappa^{2}}{2}\delta_{\Delta\Delta}u + \delta_{t-}u] \nonumber \\
  (\frac{1}{\rho_{s}Ah_{s}} + \frac{1}{\rho_{p}Hh^{2}_{p}})f  &=& \frac{c^{2}}{h_{s}}\delta_{x+}w_{0} + \frac{2}{k}\delta_{t-}w_{0} - J^{T}[-\kappa^{2}\delta_{\Delta\Delta}u + \frac{2}{k}\delta_{t-}u] \nonumber \\
f &=& (\mathcal{M})(\frac{c^{2}}{h_{s}}\delta_{x+}w_{0} + \frac{2}{k}\delta_{t-}w_{0} - ...\nonumber \\ && J^{T}[-\kappa^{2} \delta_{\Delta^{2}}u  + \frac{2}{k}\delta_{t-}u ]))\nonumber \\
\label{dvn:1DforceLoss}\end{aligned}$$

### Stiff String FD Scheme

$$\begin{aligned}
w^{n+1} &=& c^{2}k^{2} \delta_{\eta\eta}w - \kappa^{2}k^{2} \delta_{\eta\eta\eta\eta} - 2\sigma_{0}\delta_{t\cdot}u + 2\sigma_{1}\delta_{t-}\delta_{\eta\eta} + 2w - w^{n-1} + \frac{hk^2}{\rho A}Jf(n)
\label{dvn:stiffFDu}\end{aligned}$$

$$\begin{aligned}
\delta_{tt}w &= c^{2}\delta_{\eta\eta}w - \kappa^{2}\delta_{\eta\eta\eta\eta}w - 2\sigma_{0}\delta_{t\cdot}w + 2\sigma_{1}\delta_{t-}w\delta_{\eta\eta}w + 2w - w^{n-1} + \frac{hk^2}{\rho A}Jf(n) \nonumber \\
w^{n+1} &= c^{2}k^{2} \delta_{\eta\eta}w - \kappa^{2}k^{2} \delta_{\eta\eta\eta\eta} - 2k^2\sigma_{0}\delta_{t\cdot}w + 2k^2\sigma_{1}\delta_{t-}\delta_{\eta\eta}w + \frac{h}{\rho A}Jf(n) \nonumber \\
(1+k\sigma_{0})w^{n+1} &= \lambda^2 D_{\eta\eta}w - \mu^2 D_{\eta\eta\eta\eta} + 2k\sigma_{1}\delta_{\eta\eta}w + 2w - ... \nonumber \\
&(1-k\sigma_{0})w^{n-1}  - 2k\sigma_{1}\delta_{\eta\eta}w^{n-1}+ \frac{h}{\rho A}Jf(n) \nonumber \\
\underbrace{(1+k\sigma_{0})}_{A}w^{n+1} &= \underbrace{\left(\lambda^2 D_{\eta\eta} - \mu^2 D_{\eta\eta\eta\eta} + \frac{2k\sigma_{1}}{h^2}D_{\eta\eta} + 2I\right)}_{B}
- ...\nonumber \\
&\underbrace{\left((1-k\sigma_{0})w^{n-1}  + \frac{2k\sigma_{1}}{h^2}D_{\eta\eta}w^{n-1}\right)}_{C}+ \frac{hk^2}{\rho A}J_ff(n)\nonumber \\\end{aligned}$$

### Stiff String and Plate w/ Frequency Dependent Loss Coupling Force

$$\begin{aligned}
\left(\frac{1}{\rho_{s}Ah_{s}(1+k\sigma_{0s})}\bar{I} + \frac{1}{\rho_{p}Hh_{p}^{2}(1+k\sigma_{0p})}J^{t}J\right)f =& (\frac{\frac{c^2}{h}\delta_{\eta+}w_{0} - \frac{\kappa^{2}}{h_{s}^{2}}\delta_{\eta\eta}w_{1} + \frac{k}{2}w_{0}}{1+k\sigma_{0s}} - ...\\
&J^{T}\left[\frac{-\kappa^{2}\delta_{\Delta\Delta}u + 2\sigma\delta_{t-}\delta_{\Delta}u}{1+k\sigma_{0p} }\right])\end{aligned}$$

MATLAB Code {#appC}
===========

Included here is a series of MATLAB scripts from the material discussed
above. The hope is that it will provide a useful framework to begin
creating the desired instrument as well as crafting functions that will
make coding a little easier. Relatively small schemes can generate quite
large matrices; as a space saving measure it is advisable to use the
spare matrix format in MATLAB. The coefficient matrices are largely
compiled of zero values. A sparse matrix will simply declare zero cells
as empty space rather than an explicit floating point zero. This format
greatly increases computation speed and saves on memory.

Building the Laplacian and Bi-Harmonic {#appB:sec1}
--------------------------------------

The biharmonic can be a little difficult to grasp the construction of
the first few times round. It is worth noting the points that the
boundary grid points are zeroed out on lines 29-30. For this and the
laplacian the boundary points that are always have been left in. There
will be some wasted computation but the effect should be negligible
especially for larger matrices. An obvious further improvement would be
to construct the coefficient matrices without including any zero
boundary points.

### biharm.m {#cd:biharm}

### laplace.m

### fidimat.m {#cd:fidimat}

A composite function allowing for arbitrary creation of 1D and 2D finite
difference matrices

The Kirschoff Thin Plate {#cd:plateFD}
------------------------

KirschhoffPlate.m {#kirschhoffplate.m .unnumbered}
-----------------

The Coupled Plate and String {#cd:plate_string_FD}
----------------------------

CoupledPlateAndString.m {#coupledplateandstring.m .unnumbered}
-----------------------

Project Archive {#appD}
===============

This appendix contains a list of the files included in the archive, how
they are organised and operated as well as an explanation of any example
files. It is worth noting that all files were created and tested in
MATLAB r2016b running on macOS 10.11 El Capitan. Care has been taken to
ensure that each script should run, without alteration on any machine.
With changes in behaviour between versions of MATLAB it is difficult to
state, outside of the configuration above, that the scripts will run
without fault. A folder of older scripts that were created throughout
this project have been included for completeness. They are located in
the `Old Scripts` folder and should hopefully illustrate how the coding
has progressed.

Audio
-----

Some audio files have been provided to give an indication of the quality
of output from each script. Each audio file is the MIDI file and the
script that was used to generate the file. All audio files were created
with schemes run at 44.1 kHz unless stated otherwise in the filename.
MIDI files are also provided in the MIDI folder.

File Structure {#appD:sec1}
--------------

Files relating to an FDTD scheme have the same structure. The first
section refers to an ‘instrument file’, a file where all parameters
relating to the string and the plate, are set. These instrument files
also contain any flags relating to how the scheme is run, such as
whether it is plotted, analysis is run, the output style and boundary
conditions. Each of these should be clearly commented. There should be
no reason to change any of the code contained within the main files.
Each script has been set to run generating a mono or stereo audio vector
y. Scripts generating output from MIDI data reference midi files
contained within the folder `/midi`. All MIDI based scripts and the
`stringVerb.m` script produce a sound file containing the script named
and the file used to generate audio. A progress has been implemented in
the more process intensive scripts to give some indication of time left
to compute.

To run any script the basic flow is:-

$$\text{Set Parameters in Instrument File}\longrightarrow\text{Run Main Script}$$

By default the MIDI based script have been set to run using the
`chord.mid` file. This is a single D major chord and provides the best
example of model output quality for the shortest computation time. All
files that are runnable are on the top level of the folder as are the
relevant functions that are required.

Naming Convention
-----------------

All files follow roughly the same naming convention as well as logical
structure. Over the course of the project, naming of some sections has
changed however, along with the code comments, it should be clear where
these changes have taken place.

Variables relating to the string(s) are prefixed ‘`st_`’ and plate
variables with ‘`pl_`’. Other changes that have occurred are particular
to the coupling force vector which was named as `cpl_v` and `cpl_vc` for
current and previous time steps. With the introduction of Interpolation
in later code, this has simply been changed to `F` and `F1`
respectively.

Force
-----

For the code that generates a force vector through MIDI files, a
variable force system was implemented. Using [@ChaigneStringP2] as a
guide. The velocity data was scaled to between 1 and 40 newtons using a
simple linear scaling function.

$$f(x) = \frac{(d-c)(x - a)}{b-a} + c, \quad a\leq x\leq b, \quad c\leq f(x)\leq d$$

The same process was taken to then create a linear relationship between
a force in newtons $1\leq F \leq 40$ and a duration in milliseconds
$1\leq t \leq 2.5$.

$$t = \frac{(2.5-1)(F - 1)}{40-1} + 1$$

This appears to have the greatest effect on the strike model (raised
cosine). The pluck (half cosine) on the other sounds best at over 4ms
excitation duration. Any variation after that isn’t really noticeable.
This is the mode of generating force in all MIDI based scripts.

Stability
---------

Stability conditions, in particular the grid spacing, for each scheme
required altering with the addition of stiffness in the string. The
addition of frequency independent loss was first altered without any
change the grid spacing. This has been for the most part unproblematic.
With the more intricate such as the Piano models,
`multiStringGrandPiano.m` and `pianoModel_Interp.m`, the potential for
instability to arise was a little more likely. Sadly, that addition of
stability came quite late in the process and all credit has to go to
[@ATorin_PhD] as covered in section \[chapter3:sec\_freq\_loss\].

Calibration {#calibration}
-----------

The results of calibrating the string model against piano recordings are
contained in the files. `GrandPiano.m` and `multiStringGrandPiano.m`.
The process is covered in section \[chapter4:sec\_calibration\]. Sadly,
the strings were still attached to a sound board with no real means of
disconnecting the two. It would be interesting to record a completely
isolated set of piano strings but this was not feasible during this
project. The values of the calibration files still off at least some
greater resemblance to a physical instrument.

File Breakdown
--------------

This section is just a brief overview of each file, its functionality
and how it fits into the greater picture of the whole project

### accrue.m {#accrue.m .unnumbered}

There is no simple MATLAB function that will accumulate the values of a
vector. When considering the number of grid points in each system and
linear indexing, this would be a handy function to have. the `accrue.m`
script does just that. It is only used in the `multiStringPianoModel`
but it is integral for that script to function.

### biharm2.m {#biharm2.m .unnumbered}

The second incarnation of a 2D FDTD biharmonic sparse matrix generator.
The code has also been included in Appendix \[cd:biharm\] as it is
referenced throughout the paper. Should provide some details on use when
typing into MATLAB command window:`help biharm2`

### fidimat.m {#fidimat.m .unnumbered}

The `fiditmat.m` function was created for more general FDTD use. It is
sadly limited to only 1D and 2D and only clamped and simply supported
boundary conditions. A help file should be available when entering
`help fidimat`

### GrandPiano.m {#grandpiano.m .unnumbered}

The `GrandPiano.m` contains the parameters for the string derived from
calibrating with piano string recordings. It is exclusively called in
the `pianoModel.m` script.

### instPiano.m {#instpiano.m .unnumbered}

This is the instrument file for `pianoModel.m` and where the parameters
for that script can be edited. Getting the right values of the MIDI data
to relate to piano strings is somewhat inelegant. I would advise away
from editing this area (line 40) though it does provide transposition
for input.

### laplace2.m {#laplace2.m .unnumbered}

`laplace2.m` functions in exactly the same as `biharm2.m` except it
creates a laplacian FDTD sparse matrix.

### multiStringGrandPiano.m {#multistringgrandpiano.m .unnumbered}

This file is the `multiStringPianoModel.m` string calibration. Given the
way that the tiering of multiple strings was achieved, the indexing for
string parameters needed to altered.

### multiStringPianoModel.m {#multistringpianomodel.m .unnumbered}

This is a multi-strung version of the `PianoModel.m`. Rather than having
one string per note it mimics the piano in having one string for the
bass notes, 2 strings for the middle range and the remaining strings
with 3 strings per note, giving a grand total 225 strings being
modelled. The strings have been added by cell arrays given the uneven
spread across the rest of the strings. This in retrospect has become
more cumbersome then helpful. Upon creating this script again it would
likely include all string variables in a single vector.This script was
creating at a time when interpolation had yet to be implemented. As such
it is limited in to a minimum size to accommodate the number of strings.
Like the other MIDI based scripts it has been presented in state to play
the `chord.mid` file. Typically 18000 grid points are calculated.
Computation will likely take some time.

### multiStringPianoVars.m {#multistringpianovars.m .unnumbered}

The variables used for the `multiStringPianoModel.m` script. follows the
same rules as other instrument files though is sadly more cluttered at
the expense of including extra strings. Extra variables to consider are
`TT` which changes the what the tension should be by the percentage
given.

### note2hz {#note2hz .unnumbered}

This function takes a cell array of notes in scientific notation and
converts it to a vector of values in Hz. This is used in the
`stringVerb.m` script as well as the `plate_and_string_vars.m`
instrument file.

### PianoModel.m {#pianomodel.m .unnumbered}

The `PianoModel.m` seeks to imitate the make-up of a grand piano. It
uses the `instPiano` instrument file to get parameters and also
generates audio from a MIDI file. It will simulate 88 strings evenly
spaced across the y-axis of a plate. The spacing is dictated by the size
of the plate in the `pl_ctr`. The depth the strings are set into the
x-axis is set by `pl_rx` The visualisation has been difficult to get
correct and is not very intuitive. The plotting creates an image of the
plate with markers indicating the coupling points of the string. Typical
run times are 3 times slower than real-time. Plate size will reduce this
but the strings tend to be the largest burden at 5000 grid points.

### plate\_analysis.m {#plate_analysis.m .unnumbered}

This provides some plots of energy and modes against the scheme that was
run. This was used fairly earl in the process to ensure the plate was
behaving as intended.

### plate\_and\_string\_vars.m {#plate_and_string_vars.m .unnumbered}

The instrument file for the `thin_plate_stiff_string_loss.m` script

### PlateStringMIDI {#platestringmidi .unnumbered}

The instrument file for the `thin_plate_stiff_string_MIDI.m` script.

### plateVerbVars.m {#plateverbvars.m .unnumbered}

The instrument file for the `thin_plate_reverb.m` script. Contains
reference to variables not used throughout the script but that were used
in earlier implementations. These can be scene in older script folder
provided in the digital archive.

### read\_midi\_file.m {#read_midi_file.m .unnumbered}

This function was written for a previous project as a means of importing
MIDI data into a variable in MATLAB. To an extent it should be treated
as a black box and should not be altered.

### reverbParams.m {#reverbparams.m .unnumbered}

This is the instrument file for the `stringVerb.m` effect. It contains
mostly the same parameters of other instrument files with the inclusion
of choosing which notes to simulate.

### reverbStringPlate.m {#reverbstringplate.m .unnumbered}

An implementation of the reverb effect discussed in
\[chapter4:sec\_stringVerb\]. Audio is fed into the plate model and read
out from the strings. Constant power panning is applied to all string
running from left to right of the `tuning` variable. Like
`multiStringPianoModel.m`, the coupling is achieved without
interpolation. As such the number of strings that can be attached is
limited to the number of grid points in the y-axis. An error check is in
place to warn when this condition arises

### thin\_plate\_loss\_xy.m {#thin_plate_loss_xy.m .unnumbered}

A lossy Kirschhoff thin plate FDTD model. This script is self contained,
with no instrument files. The primary purpose was to produce
visualisations though a vector y still reads an output point on the
plate. playback is not automatic.

### thin\_plate\_reverb.m {#thin_plate_reverb.m .unnumbered}

A plate reverb effect, feed an audio file into a plate model and reads
out at another position.

### thin\_plate\_stiff\_string\_loss.m {#thin_plate_stiff_string_loss.m .unnumbered}

A single string connected to a thin plate model. This was the early
prototype for coupling and later used to generate visualisations. Plays
output when the `play_on` flag is `true`.

### thin\_plate\_stiff\_string\_MIDI.m {#thin_plate_stiff_string_midi.m .unnumbered}

Multiple strings are connected to a plate. Unlike the piano model, only
strings that play a note are simulated.

File List
---------

    accrue.m
    Audio
    biharm2.m
    fidimat.m
    GrandPiano.m
    instPiano.m
    laplace2.m
    midi
    multiStringGrandPiano.m
    multiStringPianoModel.m
    multiStringPianoVars.m
    note2hz.m
    OldSchemes
    pianoModel.m
    plate_analysis.m
    plate_and_string_vars.m
    plate_string_analysis.m
    PlateStringMIDI.m
    plateVerbVars.m
    read_midi_file.m
    reverbParams.m
    reverbStringPlate.m
    thin_plate_loss_xy.m
    thin_plate_reverb.m
    thin_plate_stiff_string_loss.m
    thin_plate_stiff_string_MIDI.m
