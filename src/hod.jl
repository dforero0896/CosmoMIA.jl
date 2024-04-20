#=
This file will contain an HOD inspired method for assigning random particles and hopefully avoid collapse steps
    
    The first step should be to position centrals (randomly if necessary) where the position should be inspired by
    the density distribution in the adjacent cells. This amounts to a CIC intepolation of the x, y, z positions.

    Check that clustering of centrals is somewhat reproduced.

    Second, using this particles as centrals, populate satellites using an NFW profile. 
    The Rs should be computed from the SO with the given delta value


    Fit HOD?

    =#


@inline function cic_weights(ddx, ddy, ddz, ii, jj, kk)
    (((1 - ddx) + ii * (-1 + 2 * ddx)) * 
    ((1 - ddy) + jj * (-1 + 2 * ddy)) *
    ((1 - ddz) + kk * (-1 + 2 * ddz)))
end #func

