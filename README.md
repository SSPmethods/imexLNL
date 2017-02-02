# imexLNL

From the paper:
*Implicit and Implicit-Explicit Strong Stability Preserving Runge–Kutta Methods with High Linear Order*

Naming Convention is `S, Pe, Pi, Plin, K`. For example the `6s3pe2pi4plinKin.mat` is a method with 6 stages, explicit order of 3, implicit order of 2, linear order of 4, with `K = inf`.

File Content
------------

Each `.mat` file contains the following variables:

- (A, b) : Explicit Butcher Coefficients
- (At, bt): Implicit Butcher Coefficients
- pex: Explicit nonlinear order
- pim: Implicit nonlinear order
- plin: Linear order
- s: Number of stages
- r: SSP coefficient of Explicit Method

Methods presented in Table-8:
-----------------------------
- 2s2pe2pi2plinKInf.mat
- 3s2pe2pi3plinKInf.mat
- 3s3pe3pi3plinKInf.mat
- 4s2pe2pi4plinKInf.mat
- 4s3pe3pi4plinKInf.mat
- 5s2pe2pi5plinKInf.mat
- 5s3pe3pi5plinKInf.mat
- 5s4pe4pi5plinKInf.mat
- 6s2pe2pi6plinKInf.mat
- 6s3pe3pi6plinKInf.mat
- 6s4pe4pi6plinKInf.mat
- 7s2pe2pi7plinKInf.mat
- 7s3pe3pi7plinKInf.mat
- 7s4pe4pi7plinKInf.mat

