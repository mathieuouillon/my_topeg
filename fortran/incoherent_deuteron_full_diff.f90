MODULE read_file_module
    IMPLICIT NONE
    REAL(KIND = 8) :: nk = 0
    REAL(KIND = 8) :: dnk(500) = 0
    REAL(KIND = 8) :: dk(500) = 0
END MODULE

MODULE imcff1_module
    IMPLICIT NONE
    REAL(KIND = 8) :: PAR = 0
    REAL(KIND = 8) :: D = 0
    REAL(KIND = 8) :: A = 0
    REAL(KIND = 8) :: B = 0
    REAL(KIND = 8) :: DKS = 0
    REAL(KIND = 8) :: MATR = 0
    REAL(KIND = 8) :: C(4) = 0
END MODULE

MODULE gpdee_module
    IMPLICIT NONE
    REAL(KIND = 8) :: PARE = 0
    REAL(KIND = 8) :: DE = 0
    REAL(KIND = 8) :: AE = 0
    REAL(KIND = 8) :: BE = 0
    REAL(KIND = 8) :: BETA = 0
    REAL(KIND = 8) :: DKAP = 0
    REAL(KIND = 8) :: MATRE = 0
    REAL(KIND = 8) :: CE(8) = 0
END MODULE

MODULE stuffmod
    USE iso_c_binding
    implicit none
    integer(c_int), bind(c, name="count") :: count_f
END MODULE

!     Integer in = 1->pdvcs | 2->ndvcs
!     Integer igpde = 1->(H,F1) |2->(E,F2)
!     Integer itot =1->BH only | 2->cross tot
!     Double  xk = |3-mom| of the inner nucleon. Value between [0,2] GeV/c
!     Double  costhp = cosinus of the polar angle of the inner nucleon. Value between [-1,1]
!     Double  phip =  azimuthal angle of  the inner nucleon. Value between [0,2\pi] (in radians)
SUBROUTINE CROSS_BOUNDNUCLEON_DEUT_F90(xb, d2g, dq2, ebeam, eh, phideg, xk, costhp, phip, in, igpde, itot, dlambda, crosssect)
    USE read_file_module
    USE stuffmod
    IMPLICIT REAL(KIND = 8)(A-H, O-Z)
    REAL(KIND = 8) :: alu
    REAL(KIND = 8) :: tbh2
    REAL(KIND = 8) :: dintf
    REAL(KIND = 8) :: res(3)
    REAL(KIND = 8) :: pp(3)

    REAL(KIND = 8), PARAMETER :: pi = ACOS(-1.0d0)
    REAL(KIND = 8), PARAMETER :: dnor = 1.d0 / (2.d0 * pi)**3
    REAL(KIND = 8), PARAMETER :: hc = 197.3270d0                                  !hbar c in MeV fm
    REAL(KIND = 8), PARAMETER :: grtorad = 2.d0 * pi / 360.d0
    REAL(KIND = 8), PARAMETER :: alphaelm = 1.d0 / (137.04d0)
    REAL(KIND = 8), PARAMETER :: gev2nanobarn = (hc**2 * 10.d0)
    REAL(KIND = 8), PARAMETER :: dmp = 0.9383d0
    REAL(KIND = 8), PARAMETER :: dmn = 0.9396d0
    REAL(KIND = 8), PARAMETER :: dmd = (dmn + dmp - 0.00226d0)

    alu = 0
    tbh2 = 0
    dintf = 0

    res = 0
    pp = 0

    IF (in==1) THEN
        dmass = dmp
    ELSE
        dmass = dmn
    END IF

    q = SQRT(dq2)
    dwmin = 1.d0 + dmass**2
    eps = 2.d0 * dmass * xb / q
    dscom = dmd * dmd + 2.d0 * ebeam * (eh + SQRT(eh**2 - dmd**2))
    y = dq2 / (dscom - dmd**2) / xb * (dmd / dmass)

    dcosthe = (-1.d0 - Y * eps**2 / 2.d0) / (SQRT(1.d0 + EPS**2))
    dsinthe = SQRT(1.d0 - dcosthe**2)

    q1z = -q / eps * SQRT(1.d0 + eps**2)
    dnu = q / eps

    IF (count_f < 1) THEN
        count_f = 1
        WRITE(*,*) "READ THE DAV18.DATI FILE ", dk(2)
        CALL mutexlock()
        OPEN(UNIT = 30, file = '../dav18.dati')
        DO ik = 1, 441
            READ(30, *) dk(ik), cz, cz, dnk(ik)
        END DO
        CLOSE(UNIT = 30)
        CALL mutexunlock()
    ELSE

    END IF

    phinexp = pi - phideg * grtorad
    ank0 = DNKDEU_F90(xk * 1000.d0 / hc) * (1000.d0 / hc)**3 / (4.d0 * pi)

    IF (in==1) THEN
        e3 = SQRT(dmn**2 + xk**2)
    ELSE
        e3 = SQRT(dmp**2 + xk**2)
    END IF

    po = (dmd - e3)

    pp(1) = xk
    pp(2) = ACOS(costhp)
    pp(3) = phip

    res(1) = SQRT(d2g * (-1.d0 + d2g / 4.d0 / dmass / dmass))
    res(3) = phinexp
    de2 = SQRT(dmass**2 + res(1)**2)

    !cccccccccccccccccccccFormule from https://arxiv.org/pdf/1708.00835.pdf

    px = pp(1) * SIN(pp(2)) * COS(pp(3))
    py = pp(1) * SIN(pp(2)) * SIN(pp(3))
    pz = pp(1) * COS(pp(2))

    off = po * po - pp(1) * pp(1)
    doub = 2.d0 * (po * dnu - pz * q1z)
    s = off - q * q + doub

    dnu2 = (s - dmass * dmass) / 2.d0 / SQRT(s)

    e2 = SQRT(dnu2 * dnu2 + dmass * dmass)
    tote = dnu2 + e2

    dnu1 = (tote * tote - off - q * q) / 2.d0 / tote

    qu1 = SQRT(q * q + dnu1 * dnu1)

    dtmin = -q * q - 2.d0 * dnu2 * (dnu1 - qu1)
    dtmax = -q * q - 2.d0 * dnu2 * (dnu1 + qu1)

    dphipp2 = -res(3) + pp(3)
    tetan = dacos(res(2))
    senp = SIN(tetan)

    p2 = res(1)

    dapi = pp(1)**2 + res(1)**2 + q1z * q1z + 2.d0 * (q1z) * pp(1) * COS(pp(2))
    db = -2.d0 * pp(1) * res(1) * COS(pp(2))
    dc = -2.d0 * pp(1) * res(1) * SIN(pp(2)) * COS(dphipp2)
    dd = -2.d0 * res(1) * (q1z)
    de = po - de2 + dnu

    dAg = de * de - dapi
    dbb = db + dd
    dradi = dbb * dbb - dAg * dAg + dc * dc

    IF (dradi>=0..and.de>=0.) THEN
        rad = SQRT(dradi)
        cosan1 = (dAg * dbb + ABS(dc) * rad) / (dbb * dbb + dc * dc)
        cosan2 = (dAg * dbb - ABS(dc) * rad) / (dbb * dbb + dc * dc)
        test1 = (de * de - dapi - dbb * cosan1) / dc
        test2 = (de * de - dapi - dbb * cosan2) / dc
        IF (test1>=0.) THEN
            res(2) = cosan1
            unozero = 1.d0
        ELSE
            IF (test2>=0.) THEN
                res(2) = cosan2
                unozero = 1.d0
            END IF
        END IF
    ELSE
        unozero = 0.
    END IF

    p2x = res(1) * DSIN(dacos(res(2))) * DCOS(res(3))
    p2y = res(1) * DSIN(dacos(res(2))) * DSIN(res(3))
    p2z = res(1) * res(2)

    dcospp2 = (px * p2x + py * p2y + pz * p2z) / (pp(1) * res(1))
    d2gmain = (de2 - po)**2 - ((p2x - px)**2 + (p2y - py)**2 + (p2z - pz)**2)

    IF (d2gmain>=dtmax .and. d2gmain<=dtmin) THEN
        CALL C0BHMOV_F90(IN, q, xb, d2gmain, ebeam, eps, y, de2, PP, RES, DCOMOV, po, dmass)
        CALL INTERFERENCE_F90(IN, q, xb, ebeam, eps, y, de2, d2gmain, PP, RES, ank0, DINT1, po, dmass, igpde, dlambda)
    ELSE
        dcomov = 0.
        dint1 = 0.
    END IF

    djac = res(1) / (2. * pi)**3 / ABS(4.d0 * dmass)
    cosph = COS(pp(3) - res(3))
    gjac = pi * djac / ABS(res(1) * (pp(1) * (SIN(pp(2)) * res(2) / senp * cosph - COS(pp(2))) - q1z))
    corr = gjac / po

    alu = ank0 * dint1 * corr * unozero / ank0 * dcomov * corr * unozero
    tbh2 = ank0 * dcomov * corr * unozero
    dintf = ank0 * dint1 * corr * unozero

    dkinmoto = q**2 * pi / (2.d0 * dmass * xb * xb * ebeam * ebeam)
    dnumint = (dintf * dkinmoto)
    dnumbh = (tbh2 * dkinmoto)

    crosssectbh = dnumbh * alphaelm**3 * gev2nanobarn / (2.d0 * pi)**3 * xk**2
    crosssectint = dnumint * alphaelm**3 * gev2nanobarn / (2.d0 * pi)**3 * xk**2

    IF (itot==1) THEN
        crosssect = crosssectbh
    ELSE
        crosssect = crosssectbh + crosssectint
    END IF
    RETURN
END SUBROUTINE CROSS_BOUNDNUCLEON_DEUT_F90

FUNCTION DNKDEU_F90(dkk) RESULT(dnkn)
    USE read_file_module
    IMPLICIT REAL(KIND = 8) (A-H, O-Z)
    REAL(KIND = 8) :: dkn = 0

    dkn = dkk
    dnkn = 0
    CALL LAGM_F90(dnk, 1, 441, 0.d0, 0.025d0, dkn, 1, dnkn)
    RETURN
END FUNCTION

SUBROUTINE C0BHMOV_F90(IN, q, xb, d2gv, ebeam, eps, y, de2, PP, RES, cbethe0, po, dm)
    IMPLICIT REAl(kind = 8)(A-H, O-Z)
    REAL(KIND = 8), PARAMETER :: pi = ACOS(-1.d0)
    REAL(KIND = 8) :: pp(3)
    REAL(KIND = 8) :: res(3)

    p = pp(1)
    thp = pp(2)
    php = pp(3)
    dp1 = p * p
    dcosp = COS(thp)
    senp = SIN(thp)

    p2 = res(1)
    dcosthn = res(2)
    phi = res(3)

    dcosthe = (-1.D0 - y * eps**2 / 2.D0) / (SQRT(1.D0 + eps**2))
    dsinthe = SQRT(1.D0 - dcosthe**2)
    dsinthn = SQRT(1.D0 - dcosthn**2)

    DJ = EBEAM * (DE2 - PO) + EBEAM * DSINTHE * P * SENP * COS(PHP) - EBEAM * DCOSTHE * (P2 * DCOSTHN - p * dcosp)
    DK = EBEAM * DSINTHE * P2 * DSINTHN
    DELTAF = PHP - PHI
    DCOSPP2 = DCOSTHN * DCOSP + DSINTHN * SENP * COS(DELTAF)
    Q2F = D2GV * (1.D6 / (197.31d0)**2)

    CALL FORMF_F90(IN, 5, -Q2F, F1, F2, GEP, gmP)

    P2X = res(1) * SIN(dacos(res(2))) * COS(res(3))
    P2Y = res(1) * SIN(dacos(res(2))) * SIN(res(3))
    P2Z = res(1) * res(2)
    PX = pp(1) * SIN(pp(2)) * COS(pp(3))
    PY = pp(1) * SIN(pp(2)) * SIN(pp(3))
    PZ = pp(1) * COS(pp(2))

    P1P2 = PO * DE2 - (px * p2x + py * p2y + pz * p2z)

    pr1 = (D2GV - 2.D0 * (DJ - DK * COS(PHI))) / Q**2

    pr2 = (Q**2 + 2.D0 * (DJ - DK * COS(PHI))) / Q**2

    PROPMOTO = (D2GV - 2.D0 * (DJ - DK * COS(PHI))) / Q**2 * (Q**2 + 2.D0 * (DJ - DK * COS(PHI))) / Q**2

    d2g = -0.41d0

    dkfree = SQRT((4.D0 * (1.D0 - XB) * XB + EPS**2) * &
            (-1.D0 + y + (y**2 * EPS**2) / 4.D0) * (D2G + &
            (2.D0 * Q**2 * (1.D0 - XB) + Q**2 * EPS**2 - &
                    2.D0 * Q**2 * (1.D0 - XB) * DSQRT(1.D0 + EPS**2)) / &
                    (4.D0 * (1.D0 - XB) * XB + EPS**2)) * &
            (D2G + (2.D0 * Q**2 * (1.D0 - XB) + Q**2 * EPS**2 + &
                    2.D0 * Q**2 * (1.D0 - XB) * DSQRT(1.D0 + EPS**2)) / &
                    (4.D0 * (1.D0 - XB) * XB + EPS**2))) / (2.D0 * Q**2)

    djfree = (1.D0 - Y - Y * EPS**2 / 2.D0) * (1.D0 + D2G / Q**2) - &
            (1.D0 - XB) * (2.D0 - Y) * D2G / Q**2

    dkfreel = q**2 * (y * (1.d0 + eps**2)) * dkfree

    djfreel = -q**2 / 2.d0 * (djfree / (y * (1.d0 + eps**2)) + 1.d0)

    propafree = (D2G - 2.D0 * (DJfreel - DKfreel * DCOS(PHI))) / Q**2 * &
            (Q**2 + 2.D0 * (DJfreel - DKfreel * DCOS(PHI))) / Q**2

    cbethe = (1.D0 / (D2GV**2 * Q**4 * PROPMOTO)) * &
            (-32.D0 * D2GV * DM**2 * F1**2 * Q**2 + &
                    8.D0 * D2GV * DM**2 * F1 * F2 * Q**2 + &
                    14.D0 * D2GV * DM**2 * F2**2 * Q**2 + &
                    16.D0 * D2GV * F1**2 * p1p2 * Q**2 - &
                    16.D0 * D2GV * F1 * F2 * p1p2 * Q**2 - &
                    30.D0 * D2GV * F2**2 * p1p2 * Q**2 + &
                    (12.D0 * D2GV * F2**2 * p1p2**2 * Q**2) / DM**2 + &
                    (2.D0 * D2GV * F2**2 * p1p2 * (dp1 - po**2) * Q**2) / &
                            DM**2 + 8.d0 * D2GV * F1 * F2 * (-dp1 + po**2) * Q**2 + &
                    6.D0 * D2GV * F2**2 * (-dp1 + po**2) * Q**2 - &
                    8.D0 * DM**2 * F1**2 * (1.D0 - D2GV / Q**2) * Q**4 - &
                    8.D0 * D2GV * F1 * F2 * (1.D0 - D2GV / Q**2) * Q**4 - &
                    8.D0 * D2GV * F2**2 * (1.D0 - D2GV / Q**2) * Q**4 + &
                    4.D0 * DM**2 * F2**2 * (1.D0 - D2GV / Q**2) * Q**4 + &
                    8.D0 * F1**2 * p1p2 * (1.D0 - D2GV / Q**2) * Q**4 - &
                    8.D0 * F2**2 * p1p2 * (1.D0 - D2GV / Q**2) * Q**4 + &
                    (4.D0 * D2GV * F2**2 * p1p2 * (1.D0 - D2GV / Q**2) * Q**4) / DM**2 + &
                    (4.D0 * F2**2 * p1p2**2 * (1.D0 - D2GV / Q**2) * Q**4) / DM**2 - &
                    32.D0 * F1 * F2 * Q**2 * D2GV**2 - 28.D0 * F2**2 * Q**2 * D2GV**2 + &
                    (16.D0 * F2**2 * p1p2 * Q**2 * D2GV**2) / DM**2 + &
                    (F2**2 * (1.D0 - D2GV / Q**2) * Q**4 * D2GV**2) / DM**2 + &
                    (4.D0 * F2**2 * Q**2 * D2GV**3) / DM**2 - &
                    16.D0 * DM**2 * F1**2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) - &
                    16.D0 * D2GV * F1 * F2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) - &
                    16.D0 * D2GV * F2**2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) + &
                    8.D0 * DM**2 * F2**2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) + &
                    16.D0 * F1**2 * p1p2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) - &
                    16.D0 * F2**2 * p1p2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) + &
                    (8.D0 * D2GV * F2**2 * p1p2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI)))) / DM**2 + &
                    (8.D0 * F2**2 * p1p2**2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI)))) / DM**2 + &
                    (2.D0 * F2**2 * D2GV**2.D0 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI)))) / DM**2 + &
                    8.D0 * DM**2 * F1**2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) + &
                    8.D0 * D2GV * F1 * F2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) + &
                    8.D0 * D2GV * F2**2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) - &
                    4.D0 * DM**2 * F2**2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) - &
                    8.D0 * F1**2 * p1p2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) + &
                    8.D0 * F2**2 * p1p2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) - &
                    (4.D0 * D2GV * F2**2 * p1p2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI)))) / &
                            DM**2 - (4.D0 * F2**2 * p1p2**2 * &
                    (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                    ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI)))) / &
                    DM**2 - (F2**2 * D2GV**2 * &
                    (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                    ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI)))) / &
                    DM**2 + 16.D0 * DM**2 * F1**2 * &
                    (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) + &
                    16.D0 * D2GV * F1 * F2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) - &
                    4.D0 * DM**2 * F1 * F2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) + &
                    14.D0 * D2GV * F2**2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) - &
                    7.D0 * DM**2 * F2**2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) - &
                    8.D0 * F1**2 * p1p2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) + &
                    8.D0 * F1 * F2 * p1p2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) + &
                    15.D0 * F2**2 * p1p2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) - &
                    (8 * D2GV * F2**2 * p1p2 * &
                            (Q**4 + D2GV**2 + &
                                    4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                                    8.D0 * (DJ - DK * DCOS(PHI))**2)) / DM**2 - &
                    (6.D0 * F2**2 * p1p2**2 * &
                            (Q**4 + D2GV**2 + &
                                    4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                                    8.D0 * (DJ - DK * DCOS(PHI))**2)) / DM**2 + &
                    4.D0 * F1 * F2 * (dp1 - po**2) * &
                            (Q**4 + D2GV**2 + &
                                    4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                                    8.D0 * (DJ - DK * DCOS(PHI))**2) + &
                    3.D0 * F2**2 * (dp1 - po**2) * &
                            (Q**4 + D2GV**2 + &
                                    4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                                    8.D0 * (DJ - DK * DCOS(PHI))**2) - &
                    (F2**2 * p1p2 * (dp1 - po**2) * &
                            (Q**4 + D2GV**2 + &
                                    4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                                    8.D0 * (DJ - DK * DCOS(PHI))**2)) / DM**2 - &
                    (2.D0 * F2**2 * D2GV**2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2)) / DM**2 + &
                    8.D0 * D2GV * F2**2 * Q**2 * &
                            D2GV - &
                    (8.D0 * D2GV * F2**2 * p1p2 * Q**2 * &
                            D2GV) / DM**2 + &
                    2.D0 * F2**2 * (1.D0 - D2GV / Q**2) * Q**4 * &
                            D2GV - &
                    (D2GV * F2**2 * (1.D0 - D2GV / Q**2) * Q**4 * &
                            D2GV) / DM**2 - &
                    (2.D0 * F2**2 * p1p2 * (1.D0 - D2GV / Q**2) * Q**4 * &
                            D2GV) / DM**2 - &
                    (4.D0 * F2**2 * Q**2 * D2GV**2 * &
                            D2GV) / DM**2 + &
                    4.D0 * F2**2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            D2GV - &
                    (2 * D2GV * F2**2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            D2GV) / DM**2 - &
                    (4 * F2**2 * p1p2 * (DJ - DK * DCOS(PHI)) * &
                            (Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            D2GV) / DM**2 - &
                    2.D0 * F2**2 * (D2GV - 2 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            D2GV + &
                    (D2GV * F2**2 * (D2GV - 2 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            D2GV) / DM**2 + &
                    (2.D0 * F2**2 * p1p2 * (D2GV - 2 * (DJ - DK * DCOS(PHI))) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            D2GV) / DM**2 - &
                    4.D0 * F2**2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) * &
                            D2GV + &
                    (2.D0 * D2GV * F2**2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) * &
                            D2GV) / DM**2 + &
                    (4.D0 * F2**2 * p1p2 * (Q**4 + D2GV**2 + &
                            4.D0 * (-D2GV + Q**2) * (DJ - DK * DCOS(PHI)) + &
                            8.D0 * (DJ - DK * DCOS(PHI))**2) * &
                            D2GV) / DM**2 + &
                    (16.D0 * ebeam * F1**2 * &
                            (po * Q**2 + DM * &
                                    (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                                    DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            (DJ - DK * DCOS(PHI)) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * Dcos(PHP) * DSIN(THP))) / (DM * XB) + &
                    (16.D0 * ebeam * F1 * F2 * &
                            (po * Q**2 + DM * &
                                    (D2GV - 2 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                                    DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            (DJ - DK * DCOS(PHI)) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP))) / (DM * XB) + &
                    (4.D0 * ebeam * F2**2 * &
                            (po * Q**2 + DM * &
                                    (D2GV - 2 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                                    DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            (DJ - DK * DCOS(PHI)) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP))) / (DM * XB) + &
                    (4.D0 * ebeam * F2**2 * p1p2 * &
                            (po * Q**2 + DM * &
                                    (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                                    DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            (DJ - DK * DCOS(PHI)) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * Dcos(PHP) * DSIN(THP))) / (DM**3 * XB) - &
                    32.D0 * ebeam * F1 * F2 * &
                            (D2GV - DM**2 + p1p2 + &
                                    (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                            (2.D0 * DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) - &
                    24.D0 * ebeam * F2**2 * &
                            (D2GV - DM**2 + p1p2 + &
                                    (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                            (2.D0 * DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) + &
                    (8.D0 * ebeam * F2**2 * p1p2 * &
                            (D2GV - DM**2 + p1p2 + &
                                    (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                            (2.D0 * DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP))) / DM**2 - &
                    16.D0 * ebeam**2 * F1 * F2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-po + dcosthe * p * DCOS(THP) + &
                                    p * dsinthe * Dcos(PHP) * Dsin(THP))**2 - &
                    12.D0 * ebeam**2 * F2**2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-po + dcosthe * p * DCOS(THP) + &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP))**2 + &
                    (4.D0 * ebeam**2 * F2**2 * p1p2 * &
                            (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-po + dcosthe * p * DCOS(THP) + &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP))**2) / DM**2 + &
                    (4.D0 * F1 * F2 * (po * Q**2 + &
                            2.D0 * DM * (D2GV - DM**2 + p1p2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p * (DSQRT(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))) / &
                            (DM**2 * XB**2) + &
                    (3.D0 * F2**2 * (po * Q**2 + &
                            2.D0 * DM * (D2GV - DM**2 + p1p2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSQRT(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))) / &
                            (DM**2 * XB**2) - &
                    (F2**2 * p1p2 * (po * Q**2 + &
                            2.D0 * DM * (D2GV - DM**2 + p1p2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p * (DSQRT(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))) / &
                            (DM**4 * XB**2) - &
                    (4.D0 * F1**2 * (po * Q**2 + &
                            DM * (D2GV - 2 * DM**2 + 2 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSQRT(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * Dcos(PHP) * Dsin(THP))) / &
                            (DM**2 * XB**2) - &
                    (4.D0 * F1 * F2 * (po * Q**2 + &
                            DM * (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSQRT(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))) / &
                            (DM**2 * XB**2) - &
                    (F2**2 * (po * Q**2 + &
                            DM * (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSQRT(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))) / &
                            (DM**2 * XB**2) - &
                    (F2**2 * p1p2 * (po * Q**2 + &
                            DM * (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))) / &
                            (DM**4 * XB**2) - &
                    (4.D0 * F1 * F2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (Sqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))** &
                                    2) / (DM**2 * XB**2) - &
                    (3.D0 * F2**2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSqrt(1 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))** &
                                    2) / (DM**2 * XB**2) + &
                    (F2**2 * p1p2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP))** &
                                    2) / (DM**4 * XB**2) + &
                    16.D0 * ebeam * F1 * F2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                            (-(Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (po - dcosthe * p * DCOS(THP) - &
                                            p * dsinthe * Dcos(PHP) * Dsin(THP))) + &
                    12.D0 * ebeam * F2**2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                            (-(Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (po - dcosthe * p * DCOS(THP) - &
                                            p * dsinthe * DCOS(PHP) * DSIN(THP))) - &
                    (4.D0 * ebeam * F2**2 * p1p2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                            (-(Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (po - dcosthe * p * DCOS(THP) - &
                                            p * dsinthe * DCOS(PHP) * DSIN(THP)))) / DM**2 - &
                    (8.D0 * D2GV * ebeam * F2**2 * &
                            (po * Q**2 + DM * &
                                    (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                                    DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            (DJ - DK * DCOS(PHI)) * &
                            (DE2 - dcosthe * p2 * DCOSTHN - &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN)) / (DM**3 * XB) - &
                    (12.D0 * ebeam * F2**2 * p1p2 * &
                            (po * Q**2 + DM * &
                                    (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                                    DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            (DJ - DK * DCOS(PHI)) * &
                            (DE2 - dcosthe * p2 * DCOSTHN - &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN)) / (DM**3 * XB) - &
                    32.D0 * ebeam**2 * F1**2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN) - &
                    32.D0 * ebeam**2 * F1 * F2 * (D2GV - 2 * (DJ - DK * DCOS(PHI))) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN) - &
                    8.D0 * ebeam**2 * F2**2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN) - &
                    (8.D0 * ebeam**2 * F2**2 * p1p2 * &
                            (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po - dcosthe * p * DCOS(THP) - &
                                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN)) / DM**2 - &
                    16.D0 * ebeam**2 * F1 * F2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN)**2 + &
                    4.D0 * ebeam**2 * F2**2 * (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN)**2 - &
                    (8.D0 * D2GV * ebeam**2 * F2**2 * &
                            (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN)**2) / DM**2 - &
                    (12.D0 * ebeam**2 * F2**2 * p1p2 * &
                            (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-DE2 + dcosthe * p2 * DCOSTHN + &
                                    p2 * dsinthe * DCOS(PHI) * DSINTHN)**2) / DM**2 - &
                    (F2**2 * p1p2 * (po * Q**2 + &
                            2.D0 * DM * (D2GV - DM**2 + p1p2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p2 * (Sqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**4 * XB**2) + &
                    (4 * F1 * F2 * (po * Q**2 + &
                            DM * (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p2 * (Sqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**2 * XB**2) + &
                    (2 * D2GV * F2**2 * (po * Q**2 + &
                            DM * (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p2 * (Sqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**4 * XB**2) - &
                    (F2**2 * (po * Q**2 + &
                            DM * (D2GV - 2 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**2 * XB**2) + &
                    (3.D0 * F2**2 * p1p2 * (po * Q**2 + &
                            DM * (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2) * XB + &
                            DSQRT(1.D0 + EPS**2) * p * Q**2 * DCOS(THP)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**4 * XB**2) + &
                    (8.D0 * F1**2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP)) * &
                            (DE2 * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**2 * XB**2) + &
                    (8.D0 * F1 * F2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP)) * &
                            (DE2 * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**2 * XB**2) + &
                    (2.D0 * F2**2 * (Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2.D0 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP)) * &
                            (DE2 * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**2 * XB**2) + &
                    (2.D0 * F2**2 * p1p2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (po * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOS(THP) + &
                                    2 * DM * ebeam * p * dsinthe * XB * DCOS(PHP) * DSIN(THP)) * &
                            (DE2 * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)) / &
                            (DM**4 * XB**2) - &
                    (4.D0 * F1 * F2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)**2) &
                            / (DM**2 * XB**2) - &
                    (2.D0 * D2GV * F2**2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)**2) &
                            / (DM**4 * XB**2) + &
                    (F2**2 * (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2.D0 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)**2) &
                            / (DM**2 * XB**2) - &
                    (3.D0 * F2**2 * p1p2 * (Q**2 + 2 * (DJ - DK * DCOS(PHI))) * &
                            (DE2 * (Q**2 - 2 * DM * ebeam * XB) + &
                                    p2 * (DSqrt(1.D0 + EPS**2) * Q**2 + &
                                            2.D0 * dcosthe * DM * ebeam * XB) * DCOSTHN + &
                                    2.D0 * DM * ebeam * p2 * dsinthe * XB * DCOS(PHI) * DSINTHN)**2) &
                            / (DM**4 * XB**2) + &
                    32.D0 * ebeam * F1**2 * &
                            (D2GV - DM**2 + p1p2 + &
                                    (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                            (2.D0 * DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) + &
                    32.D0 * ebeam * F1 * F2 * &
                            (D2GV - DM**2 + p1p2 + &
                                    (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                            (2.D0 * DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) + &
                    8.D0 * ebeam * F2**2 * (D2GV - DM**2 + p1p2 + &
                            (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) + &
                    (8.D0 * ebeam * F2**2 * p1p2 * &
                            (D2GV - DM**2 + p1p2 + &
                                    (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                            (2.D0 * DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN))) / DM**2 - &
                    16.D0 * ebeam * F1 * F2 * &
                            (D2GV - 2.D0 * DM**2 + 2.D0 * p1p2 - Q**2 + &
                                    (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                            (DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) + &
                    4.D0 * ebeam * F2**2 * (D2GV - 2 * DM**2 + 2 * p1p2 - Q**2 + &
                            (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (DM * XB)) * (DJ - DK * DCOS(PHI)) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) - &
                    16.D0 * ebeam * F1**2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (-(Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (po - dcosthe * p * DCOS(THP) - &
                                            p * dsinthe * DCOS(PHP) * DSIN(THP))) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) - &
                    16.D0 * ebeam * F1 * F2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (-(Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (po - dcosthe * p * DCOS(THP) - &
                                            p * dsinthe * DCOS(PHP) * DSIN(THP))) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) - &
                    4.D0 * ebeam * F2**2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (-(Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (po - dcosthe * p * DCOS(THP) - &
                                            p * dsinthe * DCOS(PHP) * DSIN(THP))) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) - &
                    (4.D0 * ebeam * F2**2 * p1p2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (-(Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (po - dcosthe * p * DCOS(THP) - &
                                            p * dsinthe * DCOS(PHP) * DSIN(THP))) * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN))) / DM**2 + &
                    16.D0 * F1**2 * (D2GV - DM**2 + p1p2 + &
                            (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB)) * &
                            ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                            (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (DE2 - &
                                            p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                            )) + 16.D0 * F1 * F2 * &
                    (D2GV - DM**2 + p1p2 + &
                            (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB)) * &
                    ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                    )) + 4.D0 * F2**2 * &
                    (D2GV - DM**2 + p1p2 + &
                            (Q**2 * (po + DSQRT(1.D0 + EPS**2) * p * DCOS(THP))) / &
                                    (2.D0 * DM * XB)) * &
                    ((1.D0 - D2GV / Q**2) * Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                    )) - 16.D0 * ebeam * F1**2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                    (po - dcosthe * p * DCOS(THP) - &
                            p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                    )) - 16.D0 * ebeam * F1 * F2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                    (po - dcosthe * p * DCOS(THP) - &
                            p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                    )) - 4.D0 * ebeam * F2**2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                    (po - dcosthe * p * DCOS(THP) - &
                            p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                    )) - (4.D0 * ebeam * F2**2 * p1p2 * (1.D0 - D2GV / Q**2) * &
                    Q**2 * (po - dcosthe * p * DCOS(THP) - &
                    p * dsinthe * DCOS(PHP) * DSIN(THP)) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + &
                                            dsinthe * DCOS(PHI) * DSINTHN)))) / DM**2 + &
                    16.D0 * ebeam * F1 * F2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) * &
                            (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (DE2 - &
                                            p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                            )) - 4.D0 * ebeam * F2**2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                    (DE2 - p2 * (dcosthe * DCOSTHN + &
                            dsinthe * DCOS(PHI) * DSINTHN)) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN) &
                                    )) + (8.D0 * D2GV * ebeam * F2**2 * (1.D0 - D2GV / Q**2) * &
                    Q**2 * (DE2 - p2 * &
                    (dcosthe * DCOSTHN + dsinthe * DCOS(PHI) * DSINTHN)) * &
                    (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                            (2.D0 * DM * XB) + &
                            ebeam * (DE2 - &
                                    p2 * (dcosthe * DCOSTHN + &
                                            dsinthe * DCOS(PHI) * DSINTHN)))) / DM**2 + &
                    (12.D0 * ebeam * F2**2 * p1p2 * (1.D0 - D2GV / Q**2) * Q**2 * &
                            (DE2 - p2 * (dcosthe * DCOSTHN + &
                                    dsinthe * DCOS(PHI) * DSINTHN)) * &
                            (-(Q**2 * (DE2 + DSQRT(1.D0 + EPS**2) * p2 * DCOSTHN)) / &
                                    (2.D0 * DM * XB) + &
                                    ebeam * (DE2 - &
                                            p2 * (dcosthe * DCOSTHN + &
                                                    dsinthe * DCOS(PHI) * DSINTHN)))) / DM**2)

    IF (ABS(propmoto)<=1.d-4) THEN
        cbethe0 = 0.d0
    ELSE
        cbethe0 = cbethe
    END IF
    RETURN
END

SUBROUTINE INTERFERENCE_F90(IN, q, xb, ebeam, eps, y, de2, d2gv, PP, RES, SF, DINT1, po, dm, igpde, dlambda)
    IMPLICIT REAl(kind = 8)(A-H, O-Z)
    REAL(KIND = 8) :: pp(3)
    REAL(KIND = 8) :: res(3)

    p = pp(1)
    thp = pp(2)
    php = pp(3)

    DP1 = P * P
    P2 = RES(1)
    DCOSTHN = RES(2)
    phi = res(3)

    DCOSTHE = (-1.D0 - Y * EPS**2 / 2.D0) / (DSQRT(1.D0 + EPS**2))
    DSINTHE = dsin(dacos(DCOSTHE))

    DSINTHN = dsin(dacos(DCOSTHN))

    DJ = EBEAM * (DE2 - PO) + EBEAM * DSINTHE * P * DSIN(THP) * DCOS(PHP) - &
            EBEAM * DCOSTHE * (P2 * DCOSTHN - p * dcos(thp))

    DK = EBEAM * DSINTHE * P2 * DSINTHN

    PR1 = (Q**2 + 2.D0 * (DJ - DK * DCOS(PHI))) / Q**2
    PR2 = (D2GV - 2.D0 * (DJ - DK * DCOS(PHI))) / Q**2

    !    ccc senza delta^2
    PR2T2 = (- 2.D0 * (DJ - DK * DCOS(PHI))) / Q**2

    PROPMOTO = PR1 * PR2

    DEN = Q**6 * D2GV * PROPMOTO

    DA = DSIN(PHI) * (2.D0 * p2 * Q**2 * DSINTHN * (PO * DSINTHE * &
            DSqrt(1.D0 + EPS**2) + p * &
            DSINTHE * DCOS(THP) - p * (DCOSTHE + DSqrt(1.D0 + EPS**2)) * &
            Dcos(PHP) * &
            Dsin(THP))) / (y * EPS**2)

    DApro = (2.D0 * p2 * Q**2 * DSINTHN * (PO * DSINTHE * &
            DSqrt(1.D0 + EPS**2) + p * &
            DSINTHE * DCOS(THP) - p * (DCOSTHE + DSqrt(1.D0 + EPS**2)) * &
            Dcos(PHP) * &
            Dsin(THP))) / (y * EPS**2)

    DLEP = Q**2 * (4.D0 * PR1**2 - 12.D0 * PR1 - 4.D0 * PR2T2**2 + 12.D0 * PR2T2)

    CALL IMCFFF1_F90(IN, q, xb, d2gv, ebeam, eps, y, de2, PP, RES, DIMF, SF, po, dm, igpde)

    IF (ABS(propmoto)<=1.d-4) THEN
        dint1 = 0.d0
    ELSE
        dint1 = DA * DLEP / DEN * dimf * dlambda
    END IF

    RETURN
END

SUBROUTINE IMCFFF1_F90(IN, q, xb, d2gv, ebeam, eps, y, de2, PP, RES, DIMF, SF, po, dm, igpde)
    USE imcff1_module
    IMPLICIT REAL(kind = 8)(A-H, O-Z)
    REAL(KIND = 8), PARAMETER :: pi = ACOS(-1.d0)
    REAL(KIND = 8) :: GPDIM(3)
    REAL(KIND = 8) :: PP(3)
    REAL(KIND = 8) :: RES(3)

    Q1Z = - Q / EPS * SQRT(1.D0 + EPS**2)

    p = pp(1)
    thp = pp(2)
    php = pp(3)
    DP1 = P * P

    P2 = RES(1)
    DCOSTHN = RES(2)

    SENP = SIN(ACOS(DCOSTHN))

    P2X = res(1) * SIN(dacos(res(2))) * COS(res(3))
    P2Y = res(1) * SIN(dacos(res(2))) * SIN(res(3))
    P2Z = res(1) * res(2)
    PX = pp(1) * SIN(pp(2)) * COS(pp(3))
    PY = pp(1) * SIN(pp(2)) * SIN(pp(3))
    PZ = pp(1) * COS(pp(2))

    P1P2 = PO * DE2 - (px * p2x + py * p2y + pz * p2z)

    !    ccccccc \xi con correzioni D2G/Q^2 al numeratore
    !    ccccccc

    xiv = (q**2 + d2gv / 2.d0) / (2.D0 * (po * Q / EPS - Q1Z * P * COS(THP)) + &
            PO**2 - DP1 - DM**2 + 2.D0 * &
            (Q / EPS * DE2 - P2 * DCOSTHN * Q1Z))

    xi = (q**2 + d2g / 2.d0) / (2.D0 * (po * Q / EPS - Q1Z * P * COS(THP)) + &
            PO**2 - DP1 - DM**2 + 2.D0 * &
            (Q / EPS * DE2 - P2 * DCOSTHN * Q1Z))

    IF (po<=0.D0 .or. xiv<0.D0) THEN
        dimf = 0.D0
    ELSE
        Q2F = D2GV * (1.D6 / (197.31d0)**2)
        CALL FORMF_F90(IN, 5, -Q2F, F1, F2, GEP, gmP)

        DO ISEGN = 1, 2
            IF (ISEGN==1) THEN
                XX = XIV
            ELSE
                XX = -XIV
            ENDIF

            CALL GPDNUCLEON_F90(XX, XIV, 1, In, q, xb, d2gv, dm, GPDU)
            CALL GPDNUCLEON_F90(XX, XIV, 2, In, q, xb, d2gv, dm, GPDD)
            CALL GPDNUCLEON_F90(XX, XIV, 3, In, q, xb, d2gv, dm, GPDS)

            CALL GPDEE_F90(XX, XIV, 1, q, d2gv, dm, GPDUE)
            CALL GPDEE_F90(XX, XIV, 2, q, d2gv, dm, GPDDE)
            CALL GPDEE_F90(XX, XIV, 3, q, d2gv, dm, GPDSE)

            gpdsu = dks * gpds

            IF (igpde==1) THEN
                IF (IN==1) THEN
                    GPDT = ((4.D0 / 9.D0 * (gpdu + gpdsu) + 1.D0 / 9.D0 * (gpdd + gpdsu) + 1.D0 / 9.D0 * gpds))
                ELSE
                    GPDT = ((1.D0 / 9.D0 * (gpdu + gpdsu) + 4.D0 / 9.D0 * (gpdd + gpdsu) + 1.D0 / 9.D0 * gpds))
                END IF
                GPDIM (ISEGN) = GPDT * f1
            ELSE
                IF (in==1) THEN
                    GPDTH = (4.D0 / 9.d0 * (gpdu + gpdsu) + 1.d0 / 9.d0 * (gpdd + gpdsu) + 1.d0 / 9.d0 * gpds)
                    GPDTE = (4.D0 / 9.d0 * (gpdue + gpdse) + 1.d0 / 9.d0 * (gpdde + gpdse) + 1.d0 / 9.d0 * gpdse)
                    dkingpde = -d2gv / (4.d0 * dm * dm) + xiv * (d2gv - 2.d0 * dm**2 + 2.d0 * p1p2) / (4.d0 * dm**2)
                    GPDT = (F1 * GPDTH + F2 * GPDTE * dkingpde)
                ELSE
                    GPDTH = (1.D0 / 9.d0 * (gpdu + gpdsu) + 4.d0 / 9.d0 * (gpdd + gpdsu) + 1.d0 / 9.d0 * gpds)
                    GPDTE = (1.D0 / 9.d0 * (gpduE + gpdsE) + 4.d0 / 9.d0 * (gpddE + gpdsE) + 1.d0 / 9.d0 * gpdsE)
                    dkingpde = -d2gv / (4.d0 * dm * dm) + xiv * (d2gv - 2.d0 * dm**2 + 2.d0 * p1p2) / (4.d0 * dm**2)
                    GPDT = (F1 * GPDTH + F2 * GPDTE * dkingpde)
                END IF
                GPDIM(Isegn) = GPDT
            END IF
        ENDDO
        DIMF = (GPDIM(1) - GPDIM(2))
    END IF
    RETURN
END

SUBROUTINE gpdnucleon_F90(X, XI, I, in, q, xb, d2g, dm, GPD)
    USE imcff1_module
    IMPLICIT REAl(kind = 8)(A-H, O-Z)
    REAL(KIND = 8) :: CONT(3)
    REAL(KIND = 8) :: coe(3, 4)
    REAL(KIND = 8) :: HIJ(10)
    REAL(KIND = 8) :: CONTAt(3)
    REAL(KIND = 8) :: HIJI(10)
    REAL(KIND = 8) :: hj(4)

    dq2 = q * q
    CALL COEFF_F90(I, q)
    CALL PARAREGGE_F90(I, q, dm, in)

    IF (i==1 .OR. i==2) THEN
        IF (xi<0.0001d0) THEN
            IF (X<=0.D0) x = -x

            conta = 0.d0

            DO J = 1, 4
                hj(J) = c(j) * x**(DFLOAT(J - 1) / 2.d0)
                conta = CONTA + HJ(J)
            ENDDO

            RIS = CONTA * x**(-D) * (1.d0 - x)**3
            GPD = RIS * EXP((B + A * LOG(1.d0 / X)) * D2G)

        ELSE
            IF (X>=-1.D0.AND.X<-XI) THEN
                GPD = 0.D0
            ELSE
                X1 = (X + XI) / (1.D0 + XI)
                X2 = (X - XI) / (1.D0 - XI)
                X3 = (X - XI) / (1.D0 + XI)
                X4 = (-X + XI) / (1.D0 + XI)

                DO II = 1, 2
                    IF (II==1) THEN
                        IF (X>=XI) THEN
                            CONT(1) = 0.D0
                        ELSE
                            CONT(1) = 0.D0
                            GO TO 12
                        ENDIF
                    ELSE
                        IF (X>=(-XI) .AND. X<XI) THEN
                            CONT(2) = 0.D0
                        ELSE
                            CONT(2) = 0.D0
                            GO TO 12
                        ENDIF
                    ENDIF

                    DO J = 1, 4
                        DMIJ = 2.D0 + DFLOAT(J - 1) / 2.D0 - D - (A * D2G)
                        IF (II==1) THEN
                            DGRA = (XI**2 - X) * (X1**DMIJ - X2**DMIJ) + DMIJ * XI * (1.D0 - X) * &
                                    (X1**DMIJ + x2**DMIJ)

                            HIJ(J) = (3.D0 * DGAMMA(DMIJ - 1.D0) * DGRA * C(J)) / (2.D0 * XI**3.D0 * &
                                    DGAMMA(DMIJ + 2.D0))
                        ELSE
                            DGRA = XI**2 - X + DMIJ * XI * (1.D0 - X)
                            HIJ(J) = (3.D0 * DGAMMA(DMIJ - 1.D0) * DGRA * C(J) * X1**DMIJ) / (2.D0 * XI**3.D0 * &
                                    DGAMMA(DMIJ + 2.D0))
                        ENDIF

                        cont(ii) = cont(ii) + hij(j)
                    ENDDO

                12          ENDDO

                GPD = (CONT(1) + CONT(2)) * DEXP(B * D2G)
            ENDIF
        ENDIF
    ELSE

        xt = dabs(x)
        IF (xi<0.00001d0) THEN

            conta = 0.d0

            DO J = 1, 4

                hj(J) = c(j) * xt**(DFLOAT(J - 1) / 2.d0)

                conta = CONTA + HJ(J)

            ENDDO

            RIS = CONTA * xt**(-(D + 1.d0)) * (1.d0 - xt)**5

            GPD = RIS * DEXP((B + A * DLOG(1.d0 / Xt)) * D2G)
        ELSE
            X1 = (XT + XI) / (1.D0 + XI)
            X2 = (XT - XI) / (1.D0 - XI)
            X3 = (XT - XI) / (1.D0 + XI)
            X4 = (-XT + XI) / (1.D0 + XI)

            DO III = 1, 2

                IF (III==1) THEN
                    IF (XT>=XI) THEN
                        CONTAT(1) = 0.D0
                    ELSE
                        CONTAT(1) = 0.D0
                        GO TO 66
                    ENDIF
                ELSE
                    IF (XT<XI) THEN
                        CONTAT(2) = 0.D0
                    ELSE
                        CONTAT(2) = 0.D0
                        GO TO 66
                    ENDIF
                ENDIF
                DO J = 1, 4

                    DMIJ = 3.d0 + dfloat(J - 1) / 2.D0 - D - (A * D2G) - 1.d0

                    IF (III==1) THEN
                        DGRA = ((DMIJ**2 + 2.D0) * (XI**2 - XT)**2 - (DMIJ**2 - 1.D0) * (1.D0 - XI**2) * &
                                (XT**2 - XI**2)) * (X1**DMIJ - X2**DMIJ) + 3.D0 * DMIJ * XI * &
                                (1.D0 - XT) * (XI**2 - XT) * (X1**DMIJ + X2**DMIJ)

                        HIJI(J) = (15.D0 * DGAMMA(DMIJ - 2.D0) * DGRA * C(J)) / (2.D0 * XI**5 * &
                                DGAMMA(DMIJ + 3.D0))

                    ELSE
                        DGRA = (X1**DMIJ) * ((DMIJ**2 + 2.D0) * (XI**2 - XT)**2.D0 + 3.D0 * DMIJ * &
                                XI * (1.D0 - XT) * (XI**2 - XT) - (DMIJ**2 - 1.D0) * (1.D0 - XI**2) * &
                                (XT**2 - XI**2)) - (X4**DMIJ) * ((DMIJ**2 + 2.D0) * (XI**2 + XT)**2.D0 + &
                                3.D0 * DMIJ * XI * (1.D0 + XT) * (XI**2 + XT) - (DMIJ**2 - 1.D0) * (1.D0 - XI**2) * &
                                (XT**2 - XI**2))

                        HIJI(J) = (15.D0 * DGAMMA(DMIJ - 2.D0) * DGRA * C(J)) / (2.D0 * XI**5.D0 * &
                                DGAMMA(DMIJ + 3.D0))

                    ENDIF

                    CONTAT(III) = CONTAT(III) + HIJI(J)

                ENDDO

            66       ENDDO

            GPD = (CONTAT(1) + CONTAT(2)) * DEXP(B * D2G)
        ENDIF
        IF (x<0) GPD = -(CONTAT(1) + CONTAT(2)) * DEXP(B * D2G)

    ENDIF
    RETURN
END

SUBROUTINE COEFF_F90(I, q)
    USE imcff1_module
    IMPLICIT REAL(KIND = 8)(A-H, O-Z)
    dq2 = q**2
    DL = DLOG(DQ2 / 4.D0)

    IF (I == 1) THEN
        C(1) = 1.52D0 + 0.248 * DL
        C(2) = 2.88D0 - 0.940 * DL
        C(3) = -0.095D0 * DL
        C(4) = 0.D0
    END IF
    IF (I == 2) THEN
        C(1) = 0.76D0 + 0.248 * DL
        C(2) = 3.11D0 - 1.36 * DL
        C(3) = -3.99D0 + 1.15 * DL
        C(4) = 0.D0
    END IF
    IF (I == 3) THEN
        C(1) = 0.123D0 + 0.0003D0 * DL
        C(2) = -0.327D0 - 0.004D0 * DL
        C(3) = 0.692D0 - 0.068D0 * DL
        C(4) = -0.486D0 + 0.038D0 * DL
    END IF
    RETURN
END


SUBROUTINE PARAREGGE_F90(i, q, dm, in)
    USE imcff1_module
    IMPLICIT REAL(KIND = 8)(A-H, O-Z)
    dq2 = q**2
    IF (in==1) THEN
        DM = 0.9383d0
    ELSE
        DM = 0.9396D0
    END IF
    DL = DLOG(DQ2 / 4.D0)
    IF(I == 1 .OR. I == 2) THEN
        D = 0.48D0
        A = 0.9D0
        B = 0.D0
    END IF
    IF(I == 3) THEN
        D = 0.10D0 + 0.06 * DL
        A = 0.15D0
        B = 2.58D0 + 0.25D0 * DLOG(DM**2 / (DQ2 + DM**2))
    END IF
    RETURN

END


SUBROUTINE FORMF_F90(NN, IFORM, DQM2, F1, F2, GE, GM)
    ! CCCCCCC
    ! CCCCCCC EVALUATE ELASTIC NUCLEON FORM FACTORS
    ! CCCCCCC DIRAC F1 AND F2, SACHS GE AND GM
    ! CCCCCCC
    ! CCCCCCC DQM2 = VALUE OF THE FOUR-MOMENTUM TRANSFER (FM-2)
    ! CCCCCCC
    ! CCCCCCC NN = 2 NEUTRON
    ! CCCCCCC NN = 1  PROTON
    ! CCCCCCC
    ! CCCCCCC IFORM = 1 BLATNIK-ZOFKO
    ! CCCCCCC IFORM = 2 GALSTER
    ! CCCCCCC IFORM = 3 HOHLER
    ! CCCCCCC IFORM = 4 GARI      ***********
    ! CCCCCCC IFORM = 5 DIPOLE
    ! CCCCCCC
    IMPLICIT REAL(kind = 8) (A-H, O-Z)
    DATA DMOMG/0.06340592198D0/, DMPHI/0.03746990965D0/, &
            DMOMGP/0.02780802576D0/, DMRHO/0.06654912153D0/, &
            DMRHOP/0.02994710470D0/, DMRHOPP/0.01853868385D0/, &
            AS/0.01837706810D0/, AV/0.0009255302930D0/, &
            DMUS/0.440D0/, DMUV/2.353D0/, BS/0.03542742485D0/, &
            BV/0.04282435971D0/, DM/4.75523795D0/
    DATA AMRO/0.776D0/, AMOM/0.784D0/, AKV/3.706D0/, AKS/-0.12D0/, &
            GFRO/0.377D0/, GFOM/0.411D0/, AKRO/6.62D0/, AKOM/0.163D0/, &
            ALA1/0.795D0/, ALA2/2.27D0/, ALQCD/0.29D0/, FI/1.02D0/

    TAU = 0.25D0 * DQM2 / (DM * DM)
    GO TO (1010, 1020, 1030, 1040, 1050), IFORM
    1010 RS = 1.D0 / ((1.D0 + DMOMG * DQM2) * (1.D0 + DMPHI * DQM2) * (1.D0 + DMOMGP * DQM2))
    RV = 1.D0 / ((1.D0 + DMRHO * DQM2) * (1.D0 + DMRHOP * DQM2) * (1.D0 + DMRHOPP * DQM2))
    IF(NN==2) THEN
        GE = (0.5D0 + AS * DQM2) * RS - (0.5D0 + AV * DQM2) * RV
        GM = (DMUS + 0.5D0 * BS * DQM2) * RS - (DMUV + 0.5D0 * BV * DQM2) * RV
    ELSE
        GE = (0.5D0 + AS * DQM2) * RS + (0.5D0 + AV * DQM2) * RV
        GM = (DMUS + 0.5D0 * BS * DQM2) * RS + (DMUV + 0.5D0 * BV * DQM2) * RV
    ENDIF
    GO TO 1500
    1020 GE = 1.D0 / (1.D0 + 0.0560D0 * DQM2)**2
    IF(NN==2) THEN
        GM = -1.913D0 * GE
        GE = 1.913D0 * TAU * GE / (1.D0 + 5.6D0 * TAU)
    ELSE
        GM = 2.793D0 * GE
    ENDIF
    GO TO 1500
    1030 DQ2G = 0.19731D0 * 0.19731D0 * DQM2
    F1RO = 0.5D0 * (0.955D0 + 0.090D0 / ((1.D0 + DQ2G / 0.355D0)**2)) / &
            (1.D0 + DQ2G / 0.536D0)
    F2RO = 0.5D0 * (5.335D0 + 0.962D0 / (1.D0 + DQ2G / 0.268D0)) / &
            (1.D0 + DQ2G / 0.603D0)
    TFO = AMOM * AMOM + DQ2G
    F1O = 0.71D0 / TFO
    F2O = -0.11D0 / TFO
    TFFI = 1.02D0 * 1.02D0 + DQ2G
    F1FI = -0.64D0 / TFFI
    F2FI = 0.13D0 / TFFI
    TFOI = 1.8D0 * 1.8D0 + DQ2G
    F1OI = -0.13D0 / TFOI
    F2OI = -0.02D0 / TFOI
    TFR1 = 1.21D0 * 1.21D0 + DQ2G
    F1R1 = 0.05D0 / TFR1
    F2R1 = -1.99D0 / TFR1
    TFR2 = 2.45D0 * 2.45D0 + DQ2G
    F1R2 = -0.52D0 / TFR2
    F2R2 = 0.20D0 / TFR2
    TFR3 = 2.95D0 * 2.95D0 + DQ2G
    F1R3 = 0.28D0 / TFR3
    F2R3 = 0.19D0 / TFR3
    F1S = F1O + F1FI + F1OI
    F1V = F1RO + F1R1 + F1R2 + F1R3
    F2S = F2O + F2FI + F2OI
    F2V = F2RO + F2R1 + F2R2 + F2R3
    IF(NN==2) THEN
        F1 = F1S - F1V
        F2 = F2S - F2V
    ELSE
        F1 = F1S + F1V
        F2 = F2S + F2V
    ENDIF
    GE = F1 - TAU * F2
    GM = F1 + F2
    RETURN
    1040 DQ2G = 0.19731D0 * 0.19731D0 * DQM2
    AMR2 = AMRO * AMRO
    AMO2 = AMOM * AMOM
    A1Q = ALA1 * ALA1
    A2Q = ALA2 * ALA2
    QCDQ = ALQCD * ALQCD
    Q2D = DQ2G * DLOG((A2Q + DQ2G) / QCDQ) / DLOG(A2Q / QCDQ)
    FL1 = A1Q / (A1Q + Q2D)
    FL2 = A2Q / (A2Q + Q2D)
    F1 = FL1 * FL2
    F2 = FL1 * FL2 * FL2
    F1V = 0.5D0 * F1 * (AMR2 * GFRO / (AMR2 + DQ2G) + 1.D0 - GFRO)
    F2V = 0.5D0 * F2 * (AMR2 * GFRO * AKRO / (AMR2 + DQ2G) + AKV - GFRO * AKRO)
    F1S = 0.5D0 * F1 * (AMO2 * GFOM / (AMO2 + DQ2G) + 1.D0 - GFOM)
    F2S = 0.5D0 * F2 * (AMO2 * GFOM * AKOM / (AMO2 + DQ2G) + AKS - GFOM * AKOM)
    IF(NN==2) THEN
        F1 = F1S - F1V
        F2 = F2S - F2V
    ELSE
        F1 = F1S + F1V
        F2 = F2S + F2V
    ENDIF
    GE = F1 - TAU * F2
    GM = F1 + F2
    RETURN
    1050 GE = 1.D0 / (1.D0 + 0.0560D0 * DQM2)**2
    IF(NN==2) THEN
        GM = -1.913D0 * GE
        GE = 0.D0
    ELSE
        GM = 2.793D0 * GE
    ENDIF
    1500 F1 = (GE + TAU * GM) / (1.D0 + TAU)
    F2 = (GM - GE) / (1.D0 + TAU)
    RETURN
END

! *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-
!
!      SUBROUTINE FOR INTERPOLATING POINTS GIVEN A NUMERICAL TABLE
!
SUBROUTINE LAGM_F90(G, J0, JF, R0, H, RIN, NIN, GIN)

    !  THE FUCNTION  G IS TABULATED  WITH STEP H FROM G(J0)= G(R0) TO G(JF) WITH (JF-J0)> 4 .
    !  GIN IS THE INTERPOLATED FUNCTION IN POINTS RINCCC  WHOSE FIRST VALUE IS
    !  RIN(IN0) CORRESPONDING TO GIN(I0). RIN(IN0),RIN(IN0+NIN-1)
    !     MUST RANGE BEETWEEN R0 AND R0+H*(JF-J0)
    IMPLICIT REAL(KIND = 8) (A-H, O-Z)
    REAL(KIND = 8) :: RIN, GIN
    DIMENSION G(1)
    N = JF - J0 + 1
    IF (N - 6 < 0) THEN
        WRITE(*, *) "PROBLEMS IN INTERPOLATING SUBROUTINE"
        STOP
    ELSE
        R = R0 + H * (N - 1)
        IF(RIN<R0 .OR. RIN>R) THEN
            WRITE(*, *) "POINTS OUT OF RANGE IN LAGM"
            STOP
        END IF
        DO I = 1, NIN
            AMR = (RIN - R0) / H
            MR = INT(AMR) + 1
            MR = MIN(MR, N - 3)
            IF (N - 6 <= 0) THEN
                P = AMR - 2.D0
                G0 = G(J0)
                L = J0 + 1
            ELSE
                L = MR + J0 - 2
                P = AMR - DFLOAT(L - J0 + 1)
                G0 = G(L - 1)
            END IF
            P3 = P - 3.D0
            P2 = P - 2.D0
            P1 = P - 1.D0
            P4 = P + 1.D0
            P5 = P + 2.D0
            P23 = P3 * P5
            P12 = P4 * P2
            PP1 = P * P1
            GIN = (PP1 * (.1D0 * P12 * (P5 * G(L + 4) - P3 * G0) + .5D0 * P23 * (P2 * G(L) - &
                    P4 * G(L + 3))) + P12 * P23 * (-P1 * G(L + 1) + P * G(L + 2))) / 12.D0
        END DO
        RETURN
    ENDIF
END


subroutine GPDEE_F90(x, xi, i, q, d2gv, dm, GPD)
    USE gpdee_module
    IMPLICIT REAL(KIND = 8)(A-H, O-Z)
    REAL(KIND = 8) :: cont(3) = 0
    REAL(KIND = 8) :: HIJ(10) = 0
    REAL(KIND = 8) :: contat(3) = 0
    REAL(KIND = 8) :: HIJI(10) = 0

    dq2 = q**2

    CALL COEFFE_F90(i)
    CALL PARAREGGEE_F90(i, dm, dq2)
    !CCCCCCC II = 1 DGLAP
    !CCCCCCC II = 2 ERBL

    IF(i == 1 .OR. i == 2) THEN
        IF (xi < 0.0000001) THEN
            IF (x <= 0.) THEN
                WRITE(6, *) "ATTENZIONE X NEGATIVE IN GPD(XI=0)", x, xi
                STOP
            ENDIF
            GPD = DGAMMA(2.d0 - dE + BETA) / (DGAMMA(1.d0 - de) * DGAMMA(1.d0 + &
                    BETA)) * DKAP * x**(-De) * (1.d0 - x)**BETA
        ELSE
            IF (X>=-1. .AND. X < -XI) then
                GPD = 0.d0
            ELSE
                X1 = (x + xi) / (1.d0 + xi)
                X2 = (x - xi) / (1.d0 - xi)
                X3 = (x - xi) / (1.d0 + xi)
                X4 = (-x + xi) / (1.d0 + xi)

                DO II = 1, 2
                    IF(II==1) THEN
                        IF (x >= xi) THEN
                            cont(1) = 0.d0
                        ELSE
                            cont(1) = 0.d0
                            GO TO 12
                        ENDIF
                    ELSE
                        IF(x >= (-xi) .AND. x < xi) THEN
                            cont(2) = 0.d0
                        ELSE
                            cont(2) = 0.d0
                            GO TO 12
                        ENDIF
                    ENDIF

                    DO J = 1, 8
                        DMIJ = 2.d0 + DFLOAT(j - 1) - DE - (AE * D2GV)
                        IF (II==1) THEN
                            DGRA = (XI**2 - X) * (X1**DMIJ - X2**DMIJ) + DMIJ * XI * &
                                    (1.d0 - X) * (X1**DMIJ + x2**DMIJ)
                            HIJ(J) = (3.d0 * DGAMMA(DMIJ - 1.d0) * DGRA * CE(J)) / &
                                    (2.d0 * XI**3.d0 * DGAMMA(DMIJ + 2.d0))
                        ELSE
                            DGRA = XI**2 - X + DMIJ * XI * (1.d0 - X)
                            HIJ(J) = (3.d0 * DGAMMA(DMIJ - 1.d0) * DGRA * CE(J) * &
                                    X1**DMIJ) / (2.d0 * XI**3.d0 * DGAMMA(DMIJ + 2.d0))
                        ENDIF
                        cont(ii) = cont(ii) + hij(j)
                    ENDDO

                12          ENDDO

                GPD = (cont(1) + cont(2)) * DEXP(BE * D2GV)
            endif
        endif
    ELSE

        xt = dabs(x)
        IF (xi<0.0000001) THEN

            !c     N_s = 0.155 (see Tab. 1 di https://arxiv.org/abs/0809.4126v1; checked with PARTONS)

            delta = - (1.1d0 + 0.06d0 * dlog(dq2 / 4.d0) - 0.0027d0 * dlog(dq2 / 4.d0)**2)
            GPD = -0.155d0 * xt**delta * (1.d0 - x)**7
        ELSE
            X1 = (XT + XI) / (1.d0 + XI)
            X2 = (XT - XI) / (1.d0 - XI)
            X3 = (XT - XI) / (1.d0 + XI)
            X4 = (-XT + XI) / (1.d0 + XI)

            DO III = 1, 2

                IF (III==1) THEN
                    IF (XT>=XI) THEN
                        contat(1) = 0.d0
                    ELSE
                        contat(1) = 0.d0
                        GO TO 66
                    ENDIF
                ELSE
                    IF(XT<XI)THEN
                        contat(2) = 0.d0
                    ELSE
                        contat(2) = 0.d0
                        GO TO 66
                    ENDIF
                ENDIF
                DO J = 1, 4

                    DMIJ = 3.d0 + dfloat(j - 1) - DE - (AE * D2GV)

                    IF (III==1) THEN
                        DGRA = ((DMIJ**2 + 2.d0) * (XI**2 - XT)**2 - (DMIJ**2 - &
                                1.d0) * (1.d0 - XI**2) * (XT**2 - XI**2)) * &
                                (X1**DMIJ - X2**DMIJ) + 3.d0 * DMIJ * XI * (1.d0 - XT) * &
                                (XI**2 - XT) * (X1**DMIJ + X2**DMIJ)

                        HIJI(J) = (15.d0 * DGAMMA(DMIJ - 2.d0) * DGRA * CE(J)) / &
                                (2.d0 * XI**5 * DGAMMA(DMIJ + 3.d0))

                    ELSE

                        ! CCCCCCCHSEAJ ERBL

                        DGRA = (X1**DMIJ) * ((DMIJ**2 + 2.d0) * (XI**2 - XT)**2 &
                                + 3.d0 * DMIJ * XI * (1.d0 - XT) * (XI**2 - XT) - &
                                (DMIJ**2 - 1.d0) * (1.d0 - XI**2) * (XT**2 - XI**2)) - &
                                (X4**DMIJ) * ((DMIJ**2 + 2.d0) * (XI**2 + XT)**2 + &
                                        3.d0 * DMIJ * XI * (1.d0 + XT) * (XI**2 + XT) - &
                                        (DMIJ**2 - 1.d0) * (1.d0 - XI**2) * (XT**2 - XI**2))

                        HIJI(J) = (15.d0 * DGAMMA(DMIJ - 2.d0) * DGRA * CE(J)) / &
                                (2.d0 * XI**5.d0 * DGAMMA(DMIJ + 3.d0))

                    ENDIF

                    contat(III) = contat(III) + HIJI(J)

                ENDDO

            66       ENDDO

            GPD = (contat(1) + contat(2)) * EXP(BE * D2GV)

        ENDIF

        IF (x < 0) GPD = -(contat(1) + contat(2)) * EXP(BE * D2GV)

    ENDIF
    RETURN
END


SUBROUTINE COEFFE_F90(I)
    USE gpdee_module
    IMPLICIT REAL(kind = 8)(A-H, O-Z)
    if (I == 1) then
        CE(1) = 2.2053d0
        CE(2) = -CE(1)
        CE(3) = 0.d0
        CE(4) = 0.d0
        CE(5) = 0.d0
        CE(6) = 0.d0
        CE(7) = 0.d0
        CE(8) = 0.d0
    ELSE IF (I == 2) THEN
        CE(1) = -3.114d0
        CE(2) = 8.096d0
        CE(3) = -6.477d0
        CE(4) = 1.295d0
        CE(5) = 0.1296d0
        CE(6) = 0.0362d0
        CE(7) = 0.014516d0
        CE(8) = 0.0070504d0
    ELSE IF (I == 3) THEN
        CE(1) = -0.155d0
        CE(2) = -2.d0 * CE(1)
        CE(3) = CE(1)
        CE(4) = 0.d0
        CE(5) = 0.d0
        CE(6) = 0.d0
        CE(7) = 0.d0
        CE(8) = 0.d0
    END IF
    RETURN
END

SUBROUTINE PARAREGGEE_F90(i, dm, dq2)
    USE gpdee_module
    IMPLICIT REAL(KIND = 8)(A-H, O-Z)
    dl = DLOG(dq2 / 4.)
    IF (i == 1) THEN
        DE = 0.48d0
        AE = 0.9d0
        BE = 0.d0
        BETA = 4.d0
        DKAP = 1.67d0
    ELSE IF (i == 2) THEN
        DE = 0.48d0
        AE = 0.9d0
        BE = 0.d0
        BETA = 5.6d0
        DKAP = -2.03d0
    ELSE IF (i == 3) THEN
        DE = 1.1d0 + 0.06d0 * dl - 0.0027d0 * dl * dl
        AE = 0.15d0
        BE = 2.58d0 + 0.25d0 * LOG(dm**2 / (dq2 + dm**2))
    END IF
    RETURN
END

