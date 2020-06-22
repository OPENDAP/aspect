program opendap_convert

    !     construct matrices to expand data on a nxn degree grid
    !     to spherical harmonics
    !     Jeroen Ritsema
    !p    PK obtained from JR, 2013

    parameter(LMX = 40)
    parameter(MXLENY = (LMX + 1)**2)
    parameter(MXWORK = (LMX + 1) * 4)
    parameter(MXMAT = (MXLENY * (MXLENY + 1)) / 2)

    dimension y(MXLENY), wk1(MXWORK), wk2(MXWORK), wk3(MXWORK)
    dimension wnorm(MXLENY)
    double precision ata(MXMAT), eigv(MXLENY)
    double precision d1(MXLENY), d2(MXLENY), d3(MXLENY), d4(MXLENY)
    double precision d5(MXLENY), d6(MXLENY), d7(MXLENY), d8(MXLENY)

    character*80 afl, evcfl, xyfl, getunx, afl2

    !   ----------------------------------------------
    !               invexpandxy.f
    !   ----------------------------------------------
    dimension atd(MXLENY)
    double precision evc(MXLENY), x(MXLENY)
    double precision f1, sum, w, damp

    character*80 grdfl, ofl, grdfl2

    !   ----------------------------------------------
    !                   sphexp.f
    !   ----------------------------------------------
    parameter (MAXD = 21)
    parameter (MAXP = 2891)
    parameter (MXMAT2 = (MAXD * (MAXD + 1) / 2))
    dimension z(MAXP + 1, 21)
    dimension d(MAXP + 1)
    dimension rawvec(MXLENY)

    double precision x2(MAXD), damp2

    character*80 sphfl
    character*120 strt


    !     call chekcl('| :r:1:input x,y(,z) file (order= lon, lat!!!)'
    !    1          //'|-a:r:1:output .a file'
    !    1          //'|-evc:r:1:output .evc file'
    !    1          //'|-lmax:r:1:maximum degree and order for expansion'
    !    1          //'|')

    read(5, *) xyfl
    read(5, *) afl
    read(5, *) evcfl
    read(5, *) lmax

    !      read(5,*) grdfl
    !      read(5,*) ofl
    !      read(5,*) afl
    !      read(5,*) evcfl

    if(lmax.gt.LMX) stop 'lmx.gt.LMX'
    leny = (lmax + 1)**2

    do i = 1, MXMAT
        ata(i) = 0.
    enddo

    open(19, file = xyfl, status = 'old')
    open(21, file = afl, status = 'unknown', form = 'unformatted')
    open(22, file = evcfl, status = 'unknown', form = 'unformatted')

    !     get normalisation
    call normylm(lmax, wnorm)

    nread = 0
    10    read(19, *, end = 100) xlon, xlat
    nread = nread + 1
    if(mod(nread, 1000).eq.0) write(6, '(i8,'' points read'')') nread
    call ylm(xlat, xlon, lmax, y, wk1, wk2, wk3)
    do k = 1, leny
        y(k) = y(k) * wnorm(k)
    enddo

    ind = 0
    do jj = 1, leny
        do ii = jj, leny
            ind = ind + 1
            ata(ind) = ata(ind) + dble(y(ii) * y(jj))
        enddo
    enddo

    write(21) xlon, xlat, (y(k), k = 1, leny)
    goto 10

    100   continue

    write(6, *) 'decomposing......'
    write(22) lmax
    call ahouse2(leny, 22, ata, d1, d2, d3, d4, d5, d6, d7, d8, eigv)

    !   ------------invexpandxy---------------

    !       19 = inpm
    !       21 = inpm.a
    !       22 = inpm.evc

    !      grdfl - inpm
    !      ofl - inpm.raw
    !      afl - inpm.a
    !      evcfl -inpm.evc

    rewind 19
    rewind 21
    rewind 22

    read(5, *) ofl
    read(5, *) dum
    damp = dble(dum)

    read(22) lmax
    leny = (lmax + 1)**2

    12    read(19, *, end = 101) xlon, xlat, grdv
    read(21) xlon2, xlat2, (y(k), k = 1, leny)
    if(xlon.ne.xlon2) stop 'inconsistent .a and .xyz file'
    if(xlat.ne.xlat2) stop 'inconsistent .a and .xyz file'
    do k = 1, leny
        atd(k) = atd(k) + grdv * y(k)
    enddo
    goto 12

    101   continue

    do i = 1, leny
        read(22) eigv(i), (evc(k), k = 1, leny)

        if(eigv(i).gt.1.d-7 * eigv(1)) then
            sum = 0.
            do j = 1, leny
                sum = sum + dble(atd(j)) * evc(j)
            enddo

            f1 = 1.d0 / (eigv(i) + damp)
            w = sum * f1
            do j = 1, leny
                x(j) = x(j) + w * evc(j)
            enddo
        endif
    enddo

    !     normaliseer harmonics
    call normylm(lmax, wnorm)
    do i = 1, leny
        x(i) = x(i) * dble(wnorm(i))
    enddo

    open(24, file = ofl, status = 'unknown')
    write(24, '(i3)') lmax
    !-- Hendrik multiplied by 0.01 .... WHY?
    write(24, '(5e16.8)') (x(i) * .01, i=1, leny)
    !--   write(24,'(5e16.8)') (x(i),i=1,leny)
    !      close(24)

    !   ------------sphexp---------------

    !      rawfl - inpm.raw
    !      dep1 - $dep1
    !      dep2 - $dep2
    !      sphfl - inpm.$iz.sph
    rewind 24
    read(5, *) dep1
    read(5, *) dep2
    read(5, *) sphfl

    if (dep2.lt.dep1) stop 'STOP dep2 is less than dep1'

        id1 = nint(2891. - dep2)
        id2 = nint(2891. - dep1)

        ! --  read raw file
        write(6, *) 'model lmax =', lmax
        np = (lmax + 1)**2
        if(lmax.gt.LMX) stop 'STOP >>>> lmax.gt.LMX'
        write(6, *) 'model np =', np
        read(24, '(5e16.8)') (rawvec(i), i = 1, np)
        close(24)

        !    Calculate the spline basis functions at a regular grid
        call splhsetup()

        rcmb = 3480.
        rmoho = 6346.
        rearth = 6371.
        r = rearth - dep
        m = MAXP / 2
        mp = MAXP + 1
        !     xd=-1.+2.*(r-rcmb)/(rmoho-rcmb)
        do ip = 1, 21
            do j = -m, m
                xd = float(j) / float(m)
                z(j + m + 1, ip) = splh(ip - 1, xd)
                !       write(20+ip,*) xd,z(j+m+1,ip)
            enddo
        enddo

        do i = 1, MXMAT
            ata(i) = 0.d0
        enddo

        leny = MAXD

        do ii = 1, mp
            ind = 0
            do ip1 = 1, leny
                do ip2 = ip1, leny
                    ind = ind + 1
                    ata(ind) = ata(ind) + dble(z(ii, ip1) * z(ii, ip2))
                enddo
            enddo
        enddo

        open(101, file = 'tmp.evc', status = 'unknown', form = 'unformatted')
        write(6, *) 'decomposing......'
        call ahouse(leny, 101, ata, d1, d2, d3, d4, d5, d6, d7, d8, eigv)
        close(101)

        open(101, file = 'tmp.evc', status = 'unknown', form = 'unformatted')

        do i = id1, id2
            d(i) = 1.
        enddo
        do ip = 1, leny
            t = 0.
            do j = 1, mp
                t = t + d(j) * z(j, ip)
            enddo
            atd(ip) = dble(t)
        enddo

        do i = 1, leny
            read(101) eigv(i), (evc(k), k = 1, leny)

            if(eigv(i).gt.1.d-7 * eigv(1)) then
                sum = 0.
                do j = 1, leny
                    sum = sum + dble(atd(j)) * evc(j)
                enddo

                f1 = 1.d0 / (eigv(i) + damp2)
                w = sum * f1
                do j = 1, leny
                    x2(j) = x2(j) + w * evc(j)
                enddo
            endif
        enddo
        do i = 1, leny
            write(6, *) 'AA', i, x2(i)
        enddo

        ind = 1
        open(26, file = sphfl, status = 'unknown')
        idp1 = 4
        idp2 = 24
        write(6, *) 'lmax,idp1,idp2=', lmax, idp1, idp2
        call wsphhead(lmax, idp1, idp2, strt, lstrt)
        write(26, *) strt(1:lstrt)
        do j = idp1, idp2
            spl = x2(j - 3)
            write(6, *) 'spl= ', spl
            ind = 1
            do k = 0, lmax
                ind1 = ind + 2 * k
                write(26, '(11e12.4)')(spl * rawvec(l), l=ind, ind1)
                ind = ind1 + 1
            enddo
        enddo
        close(26)

        do ii = 1, mp
            v = 0.
            do ip = 1, leny
                v = v + z(ii, ip) * sngl(x2(ip))
            enddo
            write(103, *) ii, v
        enddo

        99    continue

end

! ---------------------------------------------------------------------

subroutine wsphhead(lmax, idp1, idp2, strt, lstrt)
    parameter(MAXL = 100)
    parameter(MAXP = 24)
    dimension lask(0:MAXL), mask(MAXP)

    character*120 str1, str2, strt

    if(lmax.gt.MAXL) stop 'wsphhead >>> lmax too big for stment'
    if(idp2.gt.MAXP) stop 'wsphhead >>> idp2 too large'
    do i = 0, MAXL
        if(i.le.lmax) then
            lask(i) = 1
        else
            lask(i) = 0
        endif
    enddo

    do i = 1, 24
        mask(i) = 0
    enddo
    do i = idp1, idp2
        mask(i) = 1
    enddo

    write(str1, 10) lmax, (lask(i), i = 0, lmax)
    10    format(10x, i4, x, 100i1)
    write(str2, 20) 24, (mask(i), i = 1, 24)
    20    format(i4, x, 24i1)

    len1 = istlen(str1)
    len2 = istlen(str2)
    if((len1 + len2).gt.120) stop 'wsphhead >>> increase dim of strt'
    strt = str1(1:len1) // str2(1:len2)
    lstrt = len1 + len2 + 1
end

! ------------------------------------------------------------------------------

DOUBLE PRECISION FUNCTION DOT(A, J, B, K, N)
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    DOUBLE PRECISION A, B
    DIMENSION A(J, N), B(K, N)
    DOT = 0.D0
    DO 1 I = 1, N
    1 DOT = DOT + A(1, I) * B(1, I)
    RETURN
END

! ------------------------------------------------------------------------------

function istlen(string)
    character*1 null
    character*(*) string
    null = char(0)
    k = len(string)
    do 10 i = 1, k
        j = k + 1 - i
        if(string(j:j).eq.' '.or.string(j:j).eq.null) goto 10
        istlen = j
        goto 99
    10 continue
    istlen = 0
    99 return
end

! ------------------------------------------------------------------------------

subroutine splhsetup()
    parameter (MXKNT = 21)
    common/splhprm/spknt(MXKNT), qq0(MXKNT, MXKNT), qq(3, MXKNT, MXKNT)
    data spknt/&
            -1.00000, -0.78631, -0.59207, -0.41550, -0.25499&
            , -0.10909, 0.02353, 0.14409, 0.25367, 0.35329&
            , 0.44384, 0.52615, 0.60097, 0.66899, 0.73081&
            , 0.78701, 0.83810, 0.88454, 0.92675, 0.96512&
            , 1.00000&
            /
    dimension qqwk(3, MXKNT)
    do i = 1, MXKNT
        do j = 1, MXKNT
            if(i.eq.j) then
                qq0(j, i) = 1.
            else
                qq0(j, i) = 0.
            endif
        enddo
    enddo
    do i = 1, MXKNT
        call rspln(1, MXKNT, spknt(1), qq0(1, i), qq(1, 1, i), qqwk(1, 1))
    enddo
    return
end

! ------------------------------------------------
function splh(ind, x)

    parameter (MXKNT = 21)
    common/splhprm/spknt(MXKNT), qq0(MXKNT, MXKNT), qq(3, MXKNT, MXKNT)
    if(x.gt.1.or.x.lt.-1) then
        splh = 0.
        goto 10
    endif

    splh = rsple(1, MXKNT, spknt(1), qq0(1, MXKNT - ind), qq(1, 1, MXKNT - ind), x)

    10    continue
    return
end

! -------------------------------------------------

FUNCTION RSPLE(I1, I2, X, Y, Q, S)
    !
    ! C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
    !
    !   RSPLE RETURNS THE VALUE OF THE FUNCTION Y(X) EVALUATED AT POINT S
    !   USING THE CUBIC SPLINE COEFFICIENTS COMPUTED BY RSPLN AND SAVED IN
    !   Q.  IF S IS OUTSIDE THE INTERVAL (X(I1),X(I2)) RSPLE EXTRAPOLATES
    !   USING THE FIRST OR LAST INTERPOLATION POLYNOMIAL.  THE ARRAYS MUST
    !   BE DIMENSIONED AT LEAST - X(I2), Y(I2), AND Q(3,I2).
    !
    !                                                     -RPB
    DIMENSION X(*), Y(*), Q(3, *)
    DATA I/1/
    II = I2 - 1
    !   GUARANTEE I WITHIN BOUNDS.
    I = MAX0(I, I1)
    I = MIN0(I, II)
    !   SEE IF X IS INCREASING OR DECREASING.
    IF(X(I2) - X(I1))1, 2, 2
    !   X IS DECREASING.  CHANGE I AS NECESSARY.
    1    IF(S - X(I))3, 3, 4
    4    I = I - 1
    IF(I - I1)11, 6, 1
    3    IF(S - X(I + 1))5, 6, 6
    5    I = I + 1
    IF(I - II)3, 6, 7
    !   X IS INCREASING.  CHANGE I AS NECESSARY.
    2    IF(S - X(I + 1))8, 8, 9
    9    I = I + 1
    IF(I - II)2, 6, 7
    8    IF(S - X(I))10, 6, 6
    10   I = I - 1
    IF(I - I1)11, 6, 8
    7    I = II
    GO TO 6
    11   I = I1
    !   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
    6    H = S - X(I)
    RSPLE = Y(I) + H * (Q(1, I) + H * (Q(2, I) + H * Q(3, I)))
    RETURN
END

! --------------------------------------------------------

SUBROUTINE RSPLN(I1, I2, X, Y, Q, F)
    !
    ! C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
    !
    !   SUBROUTINE RSPLN COMPUTES CUBIC SPLINE INTERPOLATION COEFFICIENTS
    !   FOR Y(X) BETWEEN GRID POINTS I1 AND I2 SAVING THEM IN Q.  THE
    !   INTERPOLATION IS CONTINUOUS WITH CONTINUOUS FIRST AND SECOND
    !   DERIVITIVES.  IT AGREES EXACTLY WITH Y AT GRID POINTS AND WITH THE
    !   THREE POINT FIRST DERIVITIVES AT BOTH END POINTS (I1 AND I2).
    !   X MUST BE MONOTONIC BUT IF TWO SUCCESSIVE VALUES OF X ARE EQUAL
    !   A DISCONTINUITY IS ASSUMED AND SEPERATE INTERPOLATION IS DONE ON
    !   EACH STRICTLY MONOTONIC SEGMENT.  THE ARRAYS MUST BE DIMENSIONED AT
    !   LEAST - X(I2), Y(I2), Q(3,I2), AND F(3,I2).  F IS WORKING STORAGE
    !   FOR RSPLN.
    !                                                     -RPB
    !
    DIMENSION X(*), Y(*), Q(3, *), F(3, *), YY(3)
    EQUIVALENCE (YY(1), Y0)
    DATA SMALL/1.E-5/, YY/0., 0., 0./
    J1 = I1 + 1
    Y0 = 0.
    !   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL.
    IF(I2 - I1)13, 17, 8
    8    A0 = X(J1 - 1)
    !   SEARCH FOR DISCONTINUITIES.
    DO 3 I = J1, I2
        B0 = A0
        A0 = X(I)
        IF(ABS((A0 - B0) / AMAX1(A0, B0)).LT.SMALL) GO TO 4
    3    CONTINUE
    17   J1 = J1 - 1
    J2 = I2 - 2
    GO TO 5
    4    J1 = J1 - 1
    J2 = I - 3
    !   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
    5    IF(J2 + 1 - J1)9, 10, 11
    !   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
    10   J2 = J2 + 2
    Y0 = (Y(J2) - Y(J1)) / (X(J2) - X(J1))
    DO 15 J = 1, 3
        Q(J, J1) = YY(J)
    15   Q(J, J2) = YY(J)
    GO TO 12
    !   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
    11   A0 = 0.
    H = X(J1 + 1) - X(J1)
    H2 = X(J1 + 2) - X(J1)
    Y0 = H * H2 * (H2 - H)
    H = H * H
    H2 = H2 * H2
    !   CALCULATE DERIVITIVE AT NEAR END.
    B0 = (Y(J1) * (H - H2) + Y(J1 + 1) * H2 - Y(J1 + 2) * H) / Y0
    B1 = B0
    !   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
    DO 1 I = J1, J2
        H = X(I + 1) - X(I)
        Y0 = Y(I + 1) - Y(I)
        H2 = H * H
        HA = H - A0
        H2A = H - 2. * A0
        H3A = 2. * H - 3. * A0
        H2B = H2 * B0
        Q(1, I) = H2 / HA
        Q(2, I) = -HA / (H2A * H2)
        Q(3, I) = -H * H2A / H3A
        F(1, I) = (Y0 - H * B0) / (H * HA)
        F(2, I) = (H2B - Y0 * (2. * H - A0)) / (H * H2 * H2A)
        F(3, I) = -(H2B - 3. * Y0 * HA) / (H * H3A)
        A0 = Q(3, I)
    1    B0 = F(3, I)
    !   TAKE CARE OF LAST TWO ROWS.
    I = J2 + 1
    H = X(I + 1) - X(I)
    Y0 = Y(I + 1) - Y(I)
    H2 = H * H
    HA = H - A0
    H2A = H * HA
    H2B = H2 * B0 - Y0 * (2. * H - A0)
    Q(1, I) = H2 / HA
    F(1, I) = (Y0 - H * B0) / H2A
    HA = X(J2) - X(I + 1)
    Y0 = -H * HA * (HA + H)
    HA = HA * HA
    !   CALCULATE DERIVITIVE AT FAR END.
    Y0 = (Y(I + 1) * (H2 - HA) + Y(I) * HA - Y(J2) * H2) / Y0
    Q(3, I) = (Y0 * H2A + H2B) / (H * H2 * (H - 2. * A0))
    Q(2, I) = F(1, I) - Q(1, I) * Q(3, I)
    !   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
    DO 2 J = J1, J2
        K = I - 1
        Q(1, I) = F(3, K) - Q(3, K) * Q(2, I)
        Q(3, K) = F(2, K) - Q(2, K) * Q(1, I)
        Q(2, K) = F(1, K) - Q(1, K) * Q(3, K)
    2    I = K
    Q(1, I) = B1
    !   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
    9    J2 = J2 + 2
    DO 14 J = 1, 3
    14   Q(J, J2) = YY(J)
    !   SEE IF THIS DISCONTINUITY IS THE LAST.
    12   IF(J2 - I2)6, 13, 13
    !   NO.  GO BACK FOR MORE.
    6    J1 = J2 + 2
    IF(J1 - I2)8, 8, 7
    !   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
    7    DO 16 J = 1, 3
    16   Q(J, I2) = YY(J)
    !   FINI.
    13   RETURN
END

! ------------------------------------------------------------

SUBROUTINE AHOUSE(N, IU, C, Y, A, B, P, TA, TB, W, V, EV)
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    !   SCM VERSION
    DIMENSION A(*), B(*), P(*), TA(*), TB(*), W(*), Y(*), V(*)
    DIMENSION C(*), EV(*)
    EPS = 1.D-14
    UMEPS = 1.D0 - EPS
    TOL = 1.D-70
    JSKIP = 0
    KSKIP = 1
    NM1 = N - 1
    I = 1
    IDM1 = 0
    P(1) = 0.D0
    V(1) = 0.D0
    W(1) = 0.D0
    IF(N.LE.0) RETURN
    IF(N.GT.2) GO TO 4
    IF(N.EQ.2) GO TO 3
    EV(1) = C(1)
    Y(1) = 1.D0
    WRITE(IU) EV(1), Y(1)
    !END

    FILE IU
    REWIND IU
    RETURN
    3 A(1) = C(1)
    B(1) = -C(2)
    CKJ = C(3)
    IP1 = 2
    GO TO 215
    4 IP1 = I + 1
    NMI = N - I
    KJ = IDM1
    J = I
    5 JP1 = J + 1
    VJ = V(J)
    K = J
    LJ = N - J + 1
    JD = KJ + 1
    IF(KSKIP.EQ.1) GO TO 6
    PJ = P(J)
    WJ = W(J)
    6 KJ = KJ + 1
    CKJ = C(KJ)
    IF(KSKIP.EQ.1) GO TO 7
    DC = -(PJ * W(K) + WJ * P(K))
    CKJ = DC + CKJ
    C(KJ) = CKJ
    7 IF(J.GT.I) GO TO 14
    IF(K.GT.J) GO TO 8
    A(I) = CKJ
    K = K + 1
    GO TO 6
    8 Y(K) = 0.D0
    V(K) = CKJ
    K = K + 1
    IF(K.LE.N) GO TO 6
    JSKIP = 0
    SUM = DOT(V(JP1), 1, V(JP1), 1, LJ - 1)
    IF(SUM.LE.TOL) GO TO 10
    S = DSQRT(SUM)
    CSD = V(JP1)
    IF(CSD.LT.0.D0) S = -S
    V(JP1) = CSD + S
    C(JD + 1) = V(JP1)
    H = SUM + CSD * S
    B(I) = -S
    GO TO 12
    10 B(I) = 0.D0
    JSKIP = 1
    12 IDM1 = KJ
    IF(JSKIP.EQ.1.AND.KSKIP.EQ.1) GO TO 215
    J = JP1
    GO TO 5
    14 IF(JSKIP.EQ.0) GO TO 15
    K = K + 1
    IF(K.LE.N) GO TO 6
    J = JP1
    IF(J.LE.N) GO TO 5
    GO TO 215
    15 Y(K) = Y(K) + CKJ * VJ
    K = K + 1
    IF(K.LE.N) GO TO 6
    IF(J.EQ.N) GO TO 17
    Y(J) = Y(J) + DOT(C(JD + 1), 1, V(JP1), 1, LJ - 1)
    J = JP1
    GO TO 5

    17 SP = DOT(V(IP1), 1, Y(IP1), 1, NMI) / (H + H)
    DO 21 J = IP1, N
        W(J) = V(J)
    21 P(J) = (Y(J) - SP * V(J)) / H
    215 KSKIP = JSKIP
    I = IP1
    IF(I.LE.NM1) GO TO 4
    A(N) = CKJ
    B(NM1) = -B(NM1)
    B(N) = 0.D0
    U = DABS(A(1)) + DABS(B(1))
    DO 22 I = 2, N
    22 U = DMAX1(U, DABS(A(I)) + DABS(B(I)) + DABS(B(I - 1)))
    BD = U
    RBD = 1.D0 / U
    DO 23 I = 1, N
        W(I) = B(I)
        B(I) = (B(I) / U)**2
        A(I) = A(I) / U
        V(I) = 0.D0
    23 EV(I) = -1.D0
    U = 1.D0
    IK = 1
    NDIM = KJ
    1000 K = IK
    EL = EV(K)
    24 ELAM = .5D0 * (U + EL)
    DU = (4.D0 * DABS(ELAM) + RBD) * EPS
    IF(DABS(U - EL).LE.DU) GO TO 42
    IAG = 0
    Q = A(1) - ELAM
    IF(Q.GE.0.D0) IAG = IAG + 1
    DO 38 I = 2, N
        IF(Q.EQ.0.D0) X = DABS(W(I - 1) / BD) / EPS
        IF(Q.NE.0.D0) X = B(I - 1) / Q
        Q = A(I) - ELAM - X
        IF(Q.GE.0.D0) IAG = IAG + 1
    38 CONTINUE
    IF(IAG.GE.K) GO TO 39
    U = ELAM
    GO TO 24
    39 IF(IAG.EQ.K) GO TO 41
    M = K + 1
    DO 40 MM = M, IAG
    40 EV(MM) = ELAM
    41 EL = ELAM
    GO TO 24
    42 ELAM = BD * ELAM
    EV(K) = ELAM
    IF(IK.EQ.1) GO TO 44
    IF(ELAM.GE.EV(IK - 1)) EV(IK) = UMEPS * EV(IK - 1)
    44 I = IK
    II = 1
    L = N - 1
    DO 49 J = 1, N
    49 Y(J) = 1.D0
    50 DO 51 K = 1, N
        P(K) = 0.D0
        TB(K) = W(K)
    51 TA(K) = BD * A(K) - EV(I)
    J = 1
    DO 57 JP1 = 2, N
        IF(DABS(TA(J)).LT.DABS(W(J))) GO TO 53
        IF(TA(J).EQ.0.D0) TA(J) = EPS
        F = W(J) / TA(J)
        GO TO 55
        53 F = TA(J) / W(J)
        TA(J) = W(J)
        T = TA(JP1)
        TA(JP1) = TB(J)
        TB(J) = T
        P(J) = TB(JP1)
        TB(JP1) = 0.D0
        T = Y(J)
        Y(J) = Y(JP1)
        Y(JP1) = T
        55 TB(JP1) = TB(JP1) - F * P(J)
        TA(JP1) = TA(JP1) - F * TB(J)
        Y(JP1) = Y(JP1) - F * Y(J)
    57 J = JP1
    IF(TA(N).EQ.0.D0) TA(N) = EPS
    IF(TA(L).EQ.0.D0) TA(L) = EPS
    Y(N) = Y(N) / TA(N)
    Y(L) = (Y(L) - Y(N) * TB(L)) / TA(L)
    DO 62 J = 2, L
        K = N - J
        IF(TA(K).EQ.0.D0) TA(K) = EPS
    62 Y(K) = (Y(K) - Y(K + 1) * TB(K) - Y(K + 2) * P(K)) / TA(K)
    AY = DABS(Y(1))
    DO 63 J = 2, N
    63 AY = DMAX1(AY, DABS(Y(J)))
    DO 64 J = 1, N
    64 Y(J) = Y(J) / AY
    II = II + 1
    IF(II.LE.2) GO TO 50
    ID = NDIM - 2
    L = N - 2
    DO 68 J = 1, L
        ID = ID - J - 2
        M = N - J
        H = W(M - 1)
        IF(H.EQ.0.D0) GO TO 68
        JP1 = J + 1
        T = DOT(C(ID + 1), 1, Y(M), 1, JP1) / (H * C(ID + 1))
        KJ = ID
        DO 67 K = M, N
            KJ = KJ + 1
        67 Y(K) = Y(K) + T * C(KJ)
    68 CONTINUE
    XNORM = DSQRT(DOT(Y, 1, Y, 1, N))
    DO 70 J = 1, N
    70 Y(J) = Y(J) / XNORM
    WRITE(IU) EV(IK), (Y(J), J = 1, N)
    !     WRITE(6,901) IK,EV(IK)
    901 FORMAT(I6, 5X, D20.12)
    IK = IK + 1
    IF(IK.LE.N) GO TO 1000
END FILE IU
REWIND IU
RETURN

end