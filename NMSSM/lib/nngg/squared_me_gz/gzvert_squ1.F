* vert_squ1.F
* this file is part of the process {MNE1, MNE1} -> {0, MZ}
* generated by WriteSquaredME 17 May 2010 20:04

	subroutine gzvert_squ1
	implicit character (a-s,u-z)
	implicit double complex (t)

#include "vars.h"

        Cloop(1) = Cloop(1) + 
     -   (1/(288.D0*Pi**2)*(AbbSum39*MNE1*cint33(squ1)*MTR063(squ1)*
     -         (MTR215(1,1) + MTR216(1,1))*MTR266(squ1)) + 
     -      1/(144.D0*Pi**2)*
     -       (AbbSum39*MNE1*cint43(squ1)*MTR063(squ1)*
     -         (MTR215(1,1) + MTR216(1,1))*MTR266(squ1)))/
     -    (-MZ**2 + S)

	end

