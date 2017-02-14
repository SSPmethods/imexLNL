#Co-Optimized Implicit LNL Methods

This directory gives the **CO-Optimized** SSP LNL methods

These methods have explicit SSP RK paired with implicit DIRK that have large
stability regions The 10s6pInfKimplicitLinearPair.mat contains the coefficients
of the method in the paper that has ten stages, plin=6, pex=4, pim=3. The
10s4pInfKimplicitClassicPair.mat contains the coefficients of the method that
pairs with Ketcheson's SSPRK(10,4) and has pex=plin==4, and pim=3. The
3s3pInfKimplicitLinearPair.mat contains the coefficients of the SDIRK (but with
first diagonal element zero) that is A-stable and pairs with Shu-Osher
SSPRK(3,3)
