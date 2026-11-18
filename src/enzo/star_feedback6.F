#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK6 \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback6(nx, ny, nz,
     &               d, mu, te, ge, u, v, w,
     &               metal,
     &               idual, imetal, imulti_metals, imethod, 
     &               dt, dx, t, z, procnum,
     &               dunits, x1, vunits, t1, m_eject_frac, yield,
     &               distrad, diststep, distcells,
     &               npart, xstart, ystart, zstart, ibuff,
     &               xp, yp, zp, up, vp, wp,
     &               mp, tdp, tcp, metalf, type,
     &               particleID, nsn_timestep,
     &               exptime_in, mom_mult, kick_cap, feedback_log,
     &               use_tabfbk, minit, ergSNII, ergSNIa,
     &               itracksrc, metalSNII, metalSNIa,
     &               ntabZ, ntabAge,
     &               tabZ, tabAge, tabMass, tabMetal, tabEvents,
     &               stochastic, preSN, preSNmom, pSNntabZ, pSNntabAge,
     &               pSNtabZ, pSNtabAge, pSNtabMass, pSNtabMetal,
     &               pSNtabMom, useAGB, metalAGB)
c
c  RELEASES "STAR" PARTICLE MOMENTUM, MASS AND METALS
c
c  written by: Cassi Lochhaas, Claire Kopenhafer
c  date:      February 2025
c  	this is a copy of star_maker3mom.F that has been modifed to 
c       calculate the amount of momentum to inject from properties
c       of the host and surrounding cells, following Kimm & Cen (2014)
c       and Kimm et al. (2015).
c       The momentum is then injected into a region centered on the
c       host particle and defined by distrad and diststep.
c
c       Momentum injection is due to supernovae (Type II and
c       Type Ia) and stellar winds (optional).
c
c  INPUTS:
c
c    d     - density field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    metal - metallicity density field
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    z     - current redshift
c    dunits,x1,vunits,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imetal   - metallicity flag (0 - none, 1 - yes)
c    imulti_metals - flag to use multi metals zfield 1 and 2
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tdp      - dynamical time of zone in which particle created or
c               if exptime >= 0, a flag for whether the discrete 
c               explosion has occurred
c    tcp      - creation time of particle (-1 if not a star particle)
c    metalf   - star particle metal fraction
c    distrad  - feedback distribution radius in cells
c    diststep - distance in walking steps to deposit feedback
c    distcells - total number of cells over which to distribute feedback
c    npart    - particle array size specified by calling routine
c    m_eject_frac  - fraction of stellar mass ejected back to gas
c    yield    - fraction of stellar mass that is converted to metals
c    type     - particle type
c    particleID - particle ID number given by ParticleNumber
c    exptime  - delay time for discrete explosions (if set to -1,
c               continuous energy injection used)
c    mom_mult - multiplicative factor to boost or reduce momentum
c    kick_cap - 0 or 1 flag to limit velocity kick to 1000 km/s
c    nsn_timestep - lower limit of how many SNe per timestep below which
c                   to skip doing feedback (to save time for very low fractions
c                   of SNe per timestep per particle)
c
c    use_tabfdbk - flag to use tabular feedback
c    minit     - initial particle mass
c    ergSNII   - energy in ergs per SNII event
c    ergSNIa   - energy in ergs per SNIa event
c    itracksrc - track SNII, SNIa, & AGB metals separately?
c    metalSNII - field for SNII metals
c    metalSNIa - field for SNIa metals
c    ntabZ     - number of initial metallicities
c    ntabAge   - number of population ages
c    tabZ      - metal indicies for FB yield table
c    tabAge    - age indicies for FB yield table
c    tabMass   - table of SNe mass rates
c    tabMetal  - table of SNe metal rates
c    tabEvents - table of SNe event rates
c
c    preSN     - flag to use pre-SN feedback (from tables)
c    preSNmom  - flag to inject pre-SN momentum (default is just mass and metals)
c    pSNntabZ  - number of initial metallicities in table
c    pSNntabAge - number of population ages in table
c    pSNtabZ   - metal indices for pre-SN FB table
c    pSNtabAge - age indices for pre-SN FB table
c    pSNtabMass - table of pre-SN wind mass injection rates
c    pSNtabMetal - table of pre-SN wind metal injection rates
c    pSNtabMom  - table of pre-SN wind + Lbol momentum injection rates
c
c    useAGB    - flag to deposit AGB wind mass and metals (from tables)
c                into the star particle's host cell, with no momentum
c                or energy injection. Requires use_tabfbk.
c    metalAGB  - field for AGB wind metals (only updated if itracksrc)
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c
c
c-----------------------------------------------------------------------
      use tabular_feedback
      implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
#define NO_FORTRAN_DEBUG
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imetal, 
     &      imulti_metals, imethod, procnum
      R_PREC    d(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    ge(nx,ny,nz)
      R_PREC    metal(nx,ny,nz)
      R_PREC    mu(nx,ny,nz)
      R_PREC    dt, dx, z
      R_PREC    dunits, x1, vunits, t1, nsn_timestep, m_eject_frac
      P_PREC xstart, ystart, zstart, t
      P_PREC xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tdp(npart), tcp(npart), metalf(npart)
      INTG_PREC type(npart), particleID(npart)
      R_PREC 	 exptime_in,exptime
      INTG_PREC feedback_log, kick_cap
      R_PREC    mom_mult
      R_PREC    minit(npart)
      R_PREC    ergSNII, ergSNIa, age
      INTG_PREC itracksrc, ntabZ, ntabAge, use_tabfbk, stochastic
      INTG_PREC preSN, preSNmom, pSNntabZ, pSNntabAge, useAGB
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    tabZ(ntabZ), tabAge(ntabAge), metalAGB(nx,ny,nz)
      R_PREC    tabMass(4, ntabAge, ntabZ), tabMetal(4, ntabAge, ntabZ)
      R_PREC    tabEvents(2, ntabAge, ntabZ)
      R_PREC    pSNtabZ(pSNntabZ), pSNtabAge(pSNntabAge)
      R_PREC    pSNtabMass(pSNntabAge, pSNntabZ)
      R_PREC    pSNtabMetal(pSNntabAge, pSNntabZ)
      R_PREC    pSNtabMom(pSNntabAge, pSNntabZ)
      INTG_PREC distrad, diststep, distcells
c
c  Locals
c
      INTG_PREC n, ic, jc, kc, ip, jp, kp,
     &     i, j, k, cellstep, stepk, stepj, fbuff
      R_PREC mform, energy, energy_ii, energy_ia,
     &     m_eject, meject_ii, meject_ia, yield, minitial, xv1, xv2,
     &     mass_per_cell,
     &     mass_before, mass_dummy,
     &     kin_energy_before, mom_dummy, mom_inj, energy_inj,
     &     delta_ke, ke_injected, ke_after, kin_energy_inj,
     &     nsn, nsn_ii, nsn_ia, nsn_sto, nsn_ii_sto, nsn_ia_sto,
     &     mass_end, kin_energy_end, mom_start, mom_end, 
     &     kin_energy_dummy, e_cap, E_clip, e_res
      R_PREC xfcshift, yfcshift, zfcshift, xfc, yfc, zfc,
     &     xpos, ypos, zpos
      R_PREC ke_before(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC mom_per_cell(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC energy_per_cell(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC energy_per_cell1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC Z_avg, n_avg, Z_floor,
     &	  realmass, mu_cell, add_therm
      R_PREC n_cell, Z_cell, f_ent
      R_PREC te_sum, ke_sum, ge_sum, mom_sum, mass_sum,
     &     mass_term, vel_mag, delta_ge, add_ge
      R_PREC mz_eject, mzeject_ii, mzeject_ia
      R_PREC metal_per_cell, metal_ii_per_cell, metal_ia_per_cell
      R_PREC u1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC v1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC w1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC d1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC te1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC ge1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC metal1(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC nsn_grid(nx,ny,nz), nsn_ii_grid(nx,ny,nz), 
     &     nsn_ia_grid(nx,ny,nz), meject_grid(nx,ny,nz),
     &     mzeject_grid(nx,ny,nz), mzeject_ii_grid(nx,ny,nz),
     &     mzeject_ia_grid(nx,ny,nz), up_grid(nx,ny,nz),
     &     vp_grid(nx,ny,nz), wp_grid(nx,ny,nz)
      R_PREC pSNmass_grid(nx,ny,nz), pSNmetal_grid(nx,ny,nz),
     &      pSNmom_grid(nx,ny,nz), pSNup_grid(nx,ny,nz),
     &      pSNvp_grid(nx,ny,nz), pSNwp_grid(nx,ny,nz)
      R_PREC agbmass_grid(nx,ny,nz), agbmetal_grid(nx,ny,nz)
      R_PREC agb_mass_y, agb_metal_y
      R_PREC up_avg, vp_avg, wp_avg, up_avg_w, vp_avg_w, wp_avg_w
      R_PREC wind_mass, wind_metal, wind_mom, dt_pSN
      R_PREC ke_b_wind(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC mass_b_wind, kin_energy_b_wind, mom_b_wind,
     &       mass_per_cell_w, mom_per_cell_w, metal_per_cell_w,
     &       mass_dummy_w, kin_energy_dummy_w, mom_dummy_w, mass_end_w,
     &       kin_energy_end_w, mom_end_w, kin_energy_inj_w, add_therm_w,
     &       ke_after_w, delta_ke_w

      character*32 feedbacklog_particles
      character*32 feedbacklog_grid
      character*32 debug_filename
      character*8   post

c
c     Ensure we're not using ZEUS
c
      if (imethod .eq. 2_IKIND) then
         write(6,*) "Mechanical feedback (StarParticleFeedback = 64)" //
     &      " is not supported with the ZEUS hydro algorithm " //
     &      "(HydroMethod = 2). If you want to use this feedback, " //
     &      "use a different hydro method!"
         ERROR_MESSAGE
      endif
c
c     Initialize feedback and debug logs
c
      if (feedback_log .eq. 1_IKIND) then
         write(post,'(i8)') procnum
         do i=1,8
            if(post(i:i).eq.' ') post(i:i)='0'
         enddo
         feedbacklog_particles = 'feedbacklog_parts_' // post // '.txt'
         open(unit=4,file=feedbacklog_particles,status='unknown',
     &        position='append')
         feedbacklog_grid = 'feedbacklog_grid_' // post // '.txt'
         open(unit=5,file=feedbacklog_grid,status='unknown',
     &        position='append')
      endif
#ifdef FORTRAN_DEBUG
      write(post,'(i8)') procnum
      do i=1,8
         if(post(i:i).eq.' ') post(i:i)='0'
      enddo
      debug_filename = 'feedbackdebug_' // post // '.txt'
      open(unit=7,file=debug_filename,status='unknown',
     &        position='append')
#endif
c
c     Initialize nsn, meject, mzeject, avg particle velocity grids to zeros
c     If doing preSN feedback, initialize pSN mass, metals, momentum grids to zeros
c
      do k=1_IKIND, nz
         do j=1_IKIND, ny
            do i=1_IKIND, nx
               nsn_grid(i,j,k) = 0._RKIND
               meject_grid(i,j,k) = 0._RKIND
               mzeject_grid(i,j,k) = 0._RKIND
               mzeject_ii_grid(i,j,k) = 0._RKIND
               mzeject_ia_grid(i,j,k) = 0._RKIND
               up_grid(i,j,k) = 0._RKIND
               vp_grid(i,j,k) = 0._RKIND
               wp_grid(i,j,k) = 0._RKIND
               if (preSN .eq. 1_IKIND) then
                  pSNmass_grid(i,j,k) = 0._RKIND
                  pSNmetal_grid(i,j,k) = 0._RKIND
                  pSNmom_grid(i,j,k) = 0._RKIND
                  pSNup_grid(i,j,k) = 0._RKIND
                  pSNvp_grid(i,j,k) = 0._RKIND
                  pSNwp_grid(i,j,k) = 0._RKIND
               endif
               if (useAGB .eq. 1_IKIND) then
                  agbmass_grid(i,j,k) = 0._RKIND
                  agbmetal_grid(i,j,k) = 0._RKIND
               endif
            enddo
         enddo
      enddo
c
c     Define energy cap if kick_cap is turned on
c
c     Specific-energy cap: (2000 km/s)^2/2 = 2e16 erg/g. If all
c     thermal, this is T ~ 1e8 K, so c_s stays below ~2000 km/s too
      e_cap = 0._RKIND
      if (kick_cap .eq. 1_IKIND) e_cap = 0.5_RKIND*(2.e8_RKIND)**2
c
c-----------------------------------------------------------------------
c
c     Loop over particles to add their number of supernovae, mass and
c     metals ejection to the locations of their host cells on the grids
c
      do n=1_IKIND, npart
c
c        If particle is a star and has non-negative mass and creation time:
         if (tcp(n) .gt. 0._RKIND .and. mp(n) .gt. 0._RKIND .and. 
     &        type(n) .eq. 2_IKIND) then
c
c     Compute index of the cell that the star particle
c     resides in.
            ip = int((xp(n) - xstart)/dx) + 1_IKIND
            jp = int((yp(n) - ystart)/dx) + 1_IKIND
            kp = int((zp(n) - zstart)/dx) + 1_IKIND
c
c     Make sure particle is in grid, skip if not
            if (xp(n) .lt. xstart .or. xp(n) .gt. xstart+dx*nx .or.
     &          yp(n) .lt. ystart .or. yp(n) .gt. ystart+dx*ny .or.
     &          zp(n) .lt. zstart .or. zp(n) .gt. zstart+dx*nz) then
#ifdef FORTRAN_DEBUG
               write(7,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart
#endif
               goto 10
            endif
c
c           Calculate age of particle
            age = (t-tcp(n)) * t1/yr_s
c
c     If using pre-SN feedback, get the wind mass, metals, and momentum
c     from Starburst99 tables if the particle is young enough
c
            if (preSN .eq. 1_IKIND .and. age .lt. 4.e6_RKIND) then
c              Cap timestep at 4 Myr so winds are not integrated beyond that
               dt_pSN = min(dt, (4.e6_RKIND - age)*yr_s/t1)
               call pSN_mass(wind_mass, minit(n), metalf(n),
     &                          age, dt_pSN, t1, pSNtabMass, pSNtabZ,
     &                          pSNtabAge, pSNntabZ, pSNntabAge)
               call pSN_metal(wind_metal, minit(n), metalf(n),
     &                          age, dt_pSN, t1, pSNtabMetal, pSNtabZ,
     &                          pSNtabAge, pSNntabZ, pSNntabAge)
               if (preSNmom .eq. 1_IKIND) then
                  call pSN_mom(wind_mom, minit(n), metalf(n), age,
     &                          dt_pSN, t1, vunits, pSNtabMom, pSNtabZ,
     &                          pSNtabAge, pSNntabZ, pSNntabAge)
               else
                  wind_mom = 0._RKIND
               endif
            else
               wind_mass = 0._RKIND
               wind_metal = 0._RKIND
               wind_mom = 0._RKIND
            endif
c
c     Default AGB wind mass/metals to zero; overwritten below if
c     useAGB is set (this requires use_tabfbk, checked at startup)
c
            agb_mass_y = 0._RKIND
            agb_metal_y = 0._RKIND
c
c     If using tabulated yields, the ejected mass, energy, and metals come from
c     tables and depend on the particle's mass, age, and metallicity
c
            if (use_tabfbk .eq. 1_IKIND) then
c
#ifdef FORTRAN_DEBUG
               write(7,*) 'creation time, age, dt:', tcp(n)*t1/yr_s,
     &                        age, dt*t1/yr_s
               write(7,*) 'metalf:',metalf(n)
#endif
c
c     Integrate yield tables to get expected number of SNe in this timestep
c     and mass injection
c
               call sne_number(nsn, nsn_ii, nsn_ia,
     &                        minit(n), metalf(n), age, 
     &                        dt, t1, dx, x1, dunits,
     &                        tabEvents, tabZ, tabAge, ntabZ, ntabAge)
               call sne_mass(m_eject, meject_ii, meject_ia,
     &                        minit(n), metalf(n), age, dt, t1,
     &                        tabMass, tabZ, tabAge, ntabZ, ntabAge)
#ifdef FORTRAN_DEBUG
               write(7,*) 'nsn,nsn_ii,nsn_ia:',nsn,nsn_ii,nsn_ia
               write(7,*) 'meject_ii,meject_ia,m_eject:',
     &              meject_ii*dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              meject_ia*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              m_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND
#endif
c
c     Disallow a SN channel entirely when its integrated ejecta 
c     mass is zero or negative, so that momentum/energy can
c     never be injected onto the grid without any accompanying mass.
c
               if (meject_ii .le. 0._RKIND) nsn_ii = 0._RKIND
               if (meject_ia .le. 0._RKIND) nsn_ia = 0._RKIND
               nsn = nsn_ii + nsn_ia
c
c     Metal feedback. Similar to total ejected mass,
c     table has metal mass yield per initial population mass.
c
                  if (imetal .eq. 1_IKIND) then
                     call sne_II_metal(mzeject_ii,
     &                 minit(n), metalf(n), age,
     &                 dt, t1,
     &                 tabMetal, tabZ, tabAge, ntabZ, ntabAge)
                     call sne_Ia_metal(mzeject_ia,
     &                 minit(n), metalf(n), age,
     &                 dt, t1,
     &                 tabMetal, tabZ, tabAge, ntabZ, ntabAge)
                     mz_eject = mzeject_ii + mzeject_ia
#ifdef FORTRAN_DEBUG
               write(7,*) 'mzeject_ii,mzeject_ia,mz_eject:',
     &              mzeject_ii*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              mzeject_ia*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              mz_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND
#endif
                  endif
c
c     If including AGB wind feedback, get the AGB ejected mass and
c     metals from the tabulated yields
c
               if (useAGB .eq. 1_IKIND) then
                  call agb_mass(agb_mass_y,
     &                        minit(n), metalf(n), age, dt, t1,
     &                        tabMass, tabZ, tabAge, ntabZ, ntabAge)
                  if (imetal .eq. 1_IKIND) then
                     call agb_metal(agb_metal_y,
     &                        minit(n), metalf(n), age, dt, t1,
     &                        tabMetal, tabZ, tabAge, ntabZ, ntabAge)
                  else
                     agb_metal_y = 0._RKIND
                  endif
               else
                  agb_mass_y = 0._RKIND
                  agb_metal_y = 0._RKIND
               endif
c
c     If doing stochastic feedback, use the expected number of SNe in
c     this timestep to draw from a Poisson distribution to get an integer
c     number
c
               if (stochastic .eq. 1_IKIND) then
                  nsn_ii_sto = nsn_ii
                  nsn_ia_sto = nsn_ia
                  call draw_stochastic(nsn_ii_sto, nsn_ia_sto, dt)
                  nsn_sto = nsn_ia_sto + nsn_ii_sto
#ifdef FORTRAN_DEBUG
                  write(7,*) 'nsn_sto, nsn_ii_sto, nsn_ia_sto:',
     &               nsn_sto, nsn_ii_sto, nsn_ia_sto
#endif
c
c     Divide the expected mass and metals by the expected number in
c     this timestep to get how much mass and metals each integer SN event
c     is expected to inject, then multiply by number of SNe.
c     Divide this by initial mass of particle so that it is an ejection fraction
c     of particle mass, then multiply by the current mass and take a minimum
c     with StarMassEjectionFraction*current mass to ensure no more than 
c     a user-specified fraction of the particle's mass
c     can be ejected at once.
c
                  if (nsn_ii_sto .gt. 0._RKIND) then
                     meject_ii = min((meject_ii / nsn_ii * nsn_ii_sto)* 
     &                           mp(n)/minit(n), m_eject_frac*mp(n))
                  else
                     meject_ii = 0._RKIND
                  endif
                  if (nsn_ia_sto .gt. 0._RKIND) then
                     meject_ia = min((meject_ia / nsn_ia * nsn_ia_sto)*
     &                           mp(n)/minit(n), m_eject_frac*mp(n))
                  else
                     meject_ia = 0._RKIND
                  endif
                  m_eject = meject_ii + meject_ia
#ifdef FORTRAN_DEBUG
                  write(7,*) 'stochastic meject_ii,meject_ia,m_eject:',
     &              meject_ii*dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              meject_ia*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              m_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND
#endif
c
c     Metals done similarly to mass, but ensure no more than StarMetalYield 
c     of current mass ejected as metals
c
                  if (imetal .eq. 1_IKIND) then
                     if (nsn_ii_sto .gt. 0._RKIND) then
                        mzeject_ii = min((mzeject_ii / nsn_ii * 
     &                    nsn_ii_sto)*mp(n)/minit(n), yield*mp(n))
                     else
                        mzeject_ii = 0._RKIND
                     endif
                     if (nsn_ia_sto .gt. 0._RKIND) then
                        mzeject_ia = min((mzeject_ia / nsn_ia * 
     &                    nsn_ia_sto)*mp(n)/minit(n), yield*mp(n))
                     else
                        mzeject_ia = 0._RKIND
                     endif
                     mz_eject = mzeject_ii + mzeject_ia
#ifdef FORTRAN_DEBUG
               write(7,*) 'stochastic mzeject_ii,mzeject_ia,mz_eject:',
     &              mzeject_ii*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              mzeject_ia*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              mz_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND
#endif
                  endif
                  nsn_ii = nsn_ii_sto
                  nsn_ia = nsn_ia_sto
                  nsn = nsn_sto
               else ! end stochastic SNe
c
c     Make sure ejected mass is not larger than particle mass
c     in continuous case as well
c
                  if (m_eject .gt. mp(n)) then
c     Keep the same ratio of Type II to Type Ia masses and metals
                     if (imetal .eq. 1_IKIND) then
                        mzeject_ii = mzeject_ii/mz_eject*yield*mp(n)
                        mzeject_ia = mzeject_ia/mz_eject*yield*mp(n)
                        mz_eject = mzeject_ii + mzeject_ia
                     endif
                     meject_ii = meject_ii/m_eject * m_eject_frac*mp(n)
                     meject_ia = meject_ia/m_eject * m_eject_frac*mp(n)
                     m_eject = meject_ii + meject_ia
                  endif
               endif ! end mass and metal ejection determination
c
c
c     Skip feedback if fewer than StarFeedbackPerSNeTimestepLimit SNe in this timestep
c
               if (nsn .lt. nsn_timestep) then
                  if (wind_mass .lt. 1e-15_RKIND) then
                     goto 10
                  else
                     nsn = 0._RKIND
                  endif
               endif
c
c     Subtract ejected mass from particle
c
               mp(n) = mp(n) - m_eject
               minitial = minit(n)
               if (preSN .eq. 1_IKIND) then
                  if (wind_mass .gt. mp(n)) then
                     wind_mass = m_eject_frac*mp(n)
                     wind_metal = yield*mp(n)
                  endif
                  mp(n) = mp(n) - wind_mass
               endif
               if (useAGB .eq. 1_IKIND) then
                  if (agb_mass_y .gt. mp(n)) then
                     agb_mass_y = m_eject_frac*mp(n)
                     agb_metal_y = yield*mp(n)
                  endif
                  mp(n) = mp(n) - agb_mass_y
               endif
c
c     End tabulated feedback
c
c     If not using tabulated feedback, get supernova rate and mass and metals
c     ejection from Cen & Ostriker method.
c
            else
c
c     The star_maker2 algorithm creates a star particle instantaneously.
c     This feedback routine can either be used to also inject energy
c     instantaneously or to do feedback over a longer period of time
c     as done in star_maker3 modeled on Cen & Ostriker's method that
c     accounts for the true (unsimulated) longer formation time of
c     a stellar population.  The instananeous injection mode is 
c     intended for use with low mass star particles that produce a
c     handful of SN.
c
c     Do instantaneous injection
c
               if (exptime_in .ge. 0._RKIND) then
c              convert explosion time to code units
                  exptime = exptime_in * Myr_s / t1
c              if dynamical time is zero or delay time hasn't elapsed yet,
c              no feedback
                  if (tdp(n) .eq. 0._RKIND
     &                .or. abs(t-tcp(n)) .lt. exptime) goto 10
                  tdp(n) = 0._RKIND
                  minitial = mp(n)
                  mform = mp(n)
               else ! end instantaneous method
c
c     Do slow injection
c
c     Determine how much of a given star particle would have been 
c     turned into stars during this timestep.  Then calculate the mass
c     which should have formed during this timestep dt using the integral
c     form of the Cen & Ostriker formula.
c
                  xv1 = (t - tcp(n))/tdp(n)
#ifdef FORTRAN_DEBUG
                  write(7,*) 't, tcp(n), tdp(n), xv1:',
     &              t*t1/yr_s,tcp(n)*t1/yr_s,tdp(n)*t1/yr_s,xv1
#endif
c     If age of particle is too small or too large, feedback is zero
                  if ((xv1 .lt. 0._RKIND) .or. (xv1 .gt. 12._RKIND)) 
     &               goto 10
c
                  xv2 = (t + dt - tcp(n))/tdp(n)
c
c     First calculate the initial mass of the star particle 
c     in question.
                  minitial = mp(n) / 
     &               (1._RKIND - m_eject_frac*
     &               (1._RKIND - (1._RKIND + xv1)*exp(-xv1)))
c     
c     Then, calculate the amount of mass that would have formed in
c     this timestep.
c     
                  mform = minitial * ((1._RKIND + xv1)*exp(-xv1) - 
     &               (1._RKIND + xv2)*exp(-xv2))
                  mform = max(min(mform, mp(n)), 0._RKIND)
               endif ! end slow method
c
c     Skip if very little mass is formed.
c
               if (mform/distcells/d(ip,jp,kp) .lt. 1e-10_RKIND .and.
     &            wind_mass .lt. 1e-15_RKIND) goto 10
c     
c     Calculate the number of SNe in this timestep from the stars formed
c     in this timestep as 1 SN per 100 Msun formed. If a small fraction, skip.
c           
               realmass = mform * dunits * x1**3._RKIND*dx**3_IKIND
               realmass = realmass /SolarMass
               nsn = realmass / 100._RKIND
               if (nsn .lt. nsn_timestep) then
                  if (wind_mass .lt. 1e-15_RKIND) then
                     goto 10
                  else
                     nsn = 0._RKIND
                  endif
               endif
c
c     Subtract ejected mass from particle
c
               mp(n) = mp(n) - mform * m_eject_frac
               if (preSN .eq. 1_IKIND) mp(n) = mp(n) - wind_mass
c
c     Convert ejected mass fraction to an actual mass
c
               m_eject = mform * m_eject_frac
c
c     Calculate the metals to eject as the sum of new metals formed
c     in the star and the ejected amount of old metals
c
               if (imetal .eq. 1) mz_eject = mform * 
     &                   (yield * (1._RKIND-metalf(n))) + 
     &                   m_eject * metalf(n)
#ifdef FORTRAN_DEBUG
               write(7,*) 'mform,yield,metalf(n),m_eject',
     &             mform*dunits*(x1*dx)**3_IKIND/SolarMass,yield,
     &             metalf(n),m_eject*dunits*(x1*dx)**3_IKIND/SolarMass
               write(7,*) 'nsn,m_eject:',nsn,m_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND
#endif
            endif ! end Cen & Ostriker method
c
c	   Set center of feedback zone
c
            xfc = xp(n)
            yfc = yp(n)
            zfc = zp(n)
            fbuff = ibuff + distrad
c
c     Check bounds - if star particle is near grid edge
c     then shift center of feedback region
c
            if (xfc .lt. xstart+fbuff*dx .or. 
     &          xfc .gt. xstart+dx*nx-fbuff*dx .or.
     &          yfc .lt. ystart+fbuff*dx .or. 
     &          yfc .gt. ystart+dx*ny-fbuff*dx .or.
     &          zfc .lt. zstart+fbuff*dx .or. 
     &          zfc .gt. zstart+dx*nz-fbuff*dx) then
#ifdef FORTRAN_DEBUG
               write(7,*) 'warning1: star feedback zone shifted',
     &              xfc,yfc,zfc, xstart, ystart, zstart,fbuff
#endif
c               
	            xfcshift = xfc
	            yfcshift = yfc 
	            zfcshift	= zfc
c
	            xfc = max(xfc,xstart+fbuff*dx) 
	            yfc = max(yfc,ystart+fbuff*dx) 
	            zfc = max(zfc,zstart+fbuff*dx) 
c              
	            xfc = min(xfc,xstart+dx*nx-fbuff*dx) 
	            yfc = min(yfc,ystart+dx*ny-fbuff*dx) 
	            zfc = min(zfc,zstart+dx*nz-fbuff*dx)
c
	            xfcshift = xfcshift - xfc
	            yfcshift = yfcshift - yfc
	            zfcshift	= zfcshift - zfc
c
#ifdef FORTRAN_DEBUG
	            write(7,*) 'warning2: star feedback zone shifted',
     &              xfc,yfc,zfc,dx,xfcshift,yfcshift,zfcshift, 
     &		        nx, ny,nz
#endif
            endif ! end grid edge check
c
c     Compute index of the center of feedback zone
c     Note that this is different than star_feedback3mom because
c     this is center-indexed and star_feedback3mom is corner-indexed
c 
            xpos = (xfc - xstart)/dx
            ypos = (yfc - ystart)/dx
            zpos = (zfc - zstart)/dx
c
            ic = int(xpos) + 1_IKIND
            jc = int(ypos) + 1_IKIND
            kc = int(zpos) + 1_IKIND
c
c     Clamp feedback center so feedback region [ic-distrad, ic+distrad]
c     stays within the allocated array, regardless of floating-point
c     rounding in the position-to-index conversion above or parameter values.
c
            ic = max(fbuff+1_IKIND, min(nx-fbuff, ic))
            jc = max(fbuff+1_IKIND, min(ny-fbuff, jc))
            kc = max(fbuff+1_IKIND, min(nz-fbuff, kc))
c
c     Put nsn, meject, mzeject, and weighted particle velocity values on grid at this location
c     Weight average velocity for SNe by number of SNe going off in each particle
c
            nsn_grid(ic,jc,kc) = nsn_grid(ic,jc,kc) + nsn
            if (use_tabfbk .eq. 1_IKIND) then
               nsn_ii_grid(ic,jc,kc) = nsn_ii_grid(ic,jc,kc) + nsn_ii
               nsn_ia_grid(ic,jc,kc) = nsn_ia_grid(ic,jc,kc) + nsn_ia
            endif
            meject_grid(ic,jc,kc) = meject_grid(ic,jc,kc) + m_eject
            up_grid(ic,jc,kc) = up_grid(ic,jc,kc) + up(n)*nsn
            vp_grid(ic,jc,kc) = vp_grid(ic,jc,kc) + vp(n)*nsn
            wp_grid(ic,jc,kc) = wp_grid(ic,jc,kc) + wp(n)*nsn
            if (imetal .eq. 1)
     &       mzeject_grid(ic,jc,kc) = mzeject_grid(ic,jc,kc) + mz_eject
            if (use_tabfbk .eq. 1_IKIND) then
               mzeject_ii_grid(ic,jc,kc) = mzeject_ii_grid(ic,jc,kc) 
     &           + mzeject_ii
               mzeject_ia_grid(ic,jc,kc) = mzeject_ia_grid(ic,jc,kc)
     &           + mzeject_ia
            endif
c
c     If using pre-SN feedback, add values to grids. Weight average velocity
c     by amount of wind mass ejected by each particle.
c
            if (preSN .eq. 1_IKIND) then
               pSNmass_grid(ic,jc,kc) = pSNmass_grid(ic,jc,kc) + 
     &                                  wind_mass
               pSNmetal_grid(ic,jc,kc) = pSNmetal_grid(ic,jc,kc) + 
     &                                  wind_metal
               pSNmom_grid(ic,jc,kc) = pSNmom_grid(ic,jc,kc) + wind_mom
               pSNup_grid(ic,jc,kc) = pSNup_grid(ic,jc,kc) + 
     &                                up(n)*wind_mass
               pSNvp_grid(ic,jc,kc) = pSNvp_grid(ic,jc,kc) + 
     &                                vp(n)*wind_mass
               pSNwp_grid(ic,jc,kc) = pSNwp_grid(ic,jc,kc) +
     &                                wp(n)*wind_mass
            endif
c
c     If using AGB wind feedback, add ejected mass and metals to grids.
c
            if (useAGB .eq. 1_IKIND) then
               agbmass_grid(ic,jc,kc) = agbmass_grid(ic,jc,kc) +
     &                                  agb_mass_y
               agbmetal_grid(ic,jc,kc) = agbmetal_grid(ic,jc,kc) +
     &                                  agb_metal_y
            endif
c
c     Put particle feedback data in log. If not using pre-SN feedback,
c     those values will be zero in the log.
c
            if (feedback_log .eq. 1_IKIND) then
               write(4,*) t*t1/yr_s/1e6_RKIND, ic, jc, kc,
     &             particleID(n), age/1e6_RKIND, 
     &             mp(n)*(dunits*(x1*dx)**3._RKIND)/SolarMass,
     &             minitial*(dunits*(x1*dx)**3._RKIND)/SolarMass, nsn,
     &             m_eject*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass,
     &             mz_eject*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass,
     &             wind_mass*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass,
     &             wind_metal*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass,
     &             wind_mom*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass*(vunits/1e5_RKIND)
               endif
c
 10         continue
         endif
      enddo ! end loop over star particles
c
c-----------------------------------------------------------------------
c
c     Now, the nsn, meject, and mzeject grids contain the number of supernovae,
c     amount of ejected mass, and amount of ejected metals from ALL star particles
c     on this grid. If multiple star particles contribute to a given cell,
c     their contributions are summed before the feedback is actually injected.
c     The up, vp, and wp grids contain the weighted particle velocities
c     (once normalized by nsn) for all particles in that grid cell.
c     This will be used as the "particle velocity"
c     when translating the momentum into the particle frame.
c
c
c     Loop over the feedback grids, avoiding the edges.
c     The bounds ibuff+distrad+1 to n*-ibuff-distrad ensure that the
c     feedback region [ic-distrad, ic+distrad] always stays within the
c     allocated array for any distrad or ibuff value.
c
      do kc=ibuff+distrad+1_IKIND, nz-ibuff-distrad
         do jc=ibuff+distrad+1_IKIND, ny-ibuff-distrad
            do ic=ibuff+distrad+1_IKIND, nx-ibuff-distrad
c
c ====================================================================
c     BEGIN PRE-SUPERNOVA FEEDBACK
c
c     If injecting momentum from pre-SN feedback:
c
               add_therm_w = 0._RKIND
               if (preSN .eq. 1_IKIND .and. 
     &               pSNmass_grid(ic,jc,kc) .gt. 1e-15_RKIND) then
                  if (preSNmom .eq. 1_IKIND) then
c
c     Compute the kinetic energy in the affected region
c     before momentum is added.  This
c     is needed at the end of the calculation to update the
c     total energy (te) field.
c
               ke_b_wind = 0._RKIND
               do k = -distrad, +distrad
                  stepk = abs(k)
                  do j = -distrad, +distrad
                     stepj = stepk + abs(j)
                     do i = -distrad, +distrad
                        cellstep = stepj + abs(i)
                        if (cellstep .le. diststep) then
                           ke_b_wind(i+distrad+1, j+distrad+1,
     &                        k+distrad+1)=
     &                         0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                                 (u(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                  v(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                  w(ic+i ,jc+j ,kc+k)**2_IKIND)
                        endif
                     enddo
                  enddo
               enddo
#ifdef FORTRAN_DEBUG
               write(7,*) 'ke_b_w:',ke_b_wind
     &                  *dunits*(dx*x1)**3._RKIND*vunits**2._RKIND
#endif
c
c     Weight average velocity of particles in this cell by how much
c     wind mass they eject
c
               up_avg_w = pSNup_grid(ic,jc,kc)/pSNmass_grid(ic,jc,kc)
               vp_avg_w = pSNvp_grid(ic,jc,kc)/pSNmass_grid(ic,jc,kc)
               wp_avg_w = pSNwp_grid(ic,jc,kc)/pSNmass_grid(ic,jc,kc)
c     Convert velocities into momenta and shift to frame of average velocity
               call momentum_convert(u, v, w, d, up_avg_w, vp_avg_w,
     &                    wp_avg_w, nx, ny, nz, ic, jc, kc,
     &                    +1_IKIND, distrad, diststep)
c     Sum mass, energy, and momentum before
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc,
     &            mass_b_wind, kin_energy_b_wind, mom_b_wind,
     &            distrad, diststep)
#ifdef FORTRAN_DEBUG
               write(7,*) 'up_w, vp_w, wp_w',
     &            up_avg_w*vunits/1e5_IKIND,
     &            vp_avg_w*vunits/1e5_IKIND,
     &            wp_avg_w*vunits/1e5_IKIND
               write(7,*) "m before wind:", mass_b_wind*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
               write(7,*) "ke before wind:", kin_energy_b_wind*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               write(7,*) 'mom_b_wind:',mom_b_wind*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
c           mass_before has units of mass/volume
c           kin_energy_before has units of mass*velocity^2/volume
c           mom_start has units of mass*velocity/volume
c
#endif
c
c     Zero dummy velocity, energy, and metallicity fields
c
               do k = -distrad, distrad
                  do j = -distrad, distrad
                     do i = -distrad, distrad
                        u1(i+distrad+1,j+distrad+1,k+distrad+1)=0._RKIND
                        v1(i+distrad+1,j+distrad+1,k+distrad+1)=0._RKIND
                        w1(i+distrad+1,j+distrad+1,k+distrad+1)=0._RKIND
                        d1(i+distrad+1,j+distrad+1,k+distrad+1) =
     &                     d(i+ic,j+jc,k+kc)
                        ge1(i+distrad+1,j+distrad+1,k+distrad+1) = 
     &                     0._RKIND
                        te1(i+distrad+1,j+distrad+1,k+distrad+1) = 
     &                     0._RKIND
                        metal1(i+distrad+1,j+distrad+1,k+distrad+1) =
     &                     0._RKIND
                     enddo
                  enddo
               enddo
c
c     Now add mass and momentum terms to local dummy fields
c
               mass_per_cell_w = pSNmass_grid(ic,jc,kc)/distcells
c           wind momentum per cell will be normalized by mass in that cell
               mom_per_cell_w = pSNmom_grid(ic,jc,kc)
               metal_per_cell_w = pSNmetal_grid(ic,jc,kc)/distcells
               e_res = 0._RKIND
               call add_preSN_feedback(u1, v1, w1, d1, ge1, te1, metal1,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        distrad+1_IKIND, distrad+1_IKIND,
     &                        distrad+1_IKIND,
     &                        0_IKIND, imulti_metals, idual, itracksrc,
     &                        metalSNII, metalSNIa, useAGB, metalAGB,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell_w, mom_per_cell_w,
     &                        metal_per_cell_w, dx, x1, dunits, vunits,
     &                        distrad, diststep, e_res)
c
c     Calculate mass and energy in dummy fields
c
               call sum_energy_mom(u1, v1, w1, d1,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        distrad+1_IKIND, distrad+1_IKIND,
     &                        distrad+1_IKIND,
     &                        mass_dummy_w, kin_energy_dummy_w,
     &                        mom_dummy_w, distrad, diststep)
#ifdef FORTRAN_DEBUG
               write(7,*) "m dummy wind:", mass_dummy_w*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
               write(7,*) "ke dummy wind:", kin_energy_dummy_w*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               write(7,*) 'mom_dummy wind:',mom_dummy_w*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
#endif
c
c     Now deposit the momentum, mass, and metals onto the actual grid.
c
               call add_preSN_feedback(u, v, w, d, ge, te, metal,
     &                         nx, ny, nz, ic, jc, kc,
     &                        imetal, imulti_metals, idual, itracksrc,
     &                        metalSNII, metalSNIa, useAGB, metalAGB,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell_w, mom_per_cell_w,
     &                        metal_per_cell_w,
     &                        dx, x1, dunits, vunits,
     &                        distrad, diststep, e_res)
c
c     Calculate the mass and energy in the affected region again.
c
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &                         ic, jc, kc,
     &                         mass_end_w, kin_energy_end_w, mom_end_w,
     &                         distrad, diststep)
#ifdef FORTRAN_DEBUG
               write(7,*) "m end wind:", mass_end_w*
     &               (dunits*(dx*x1)**3._RKIND)/SolarMass
               write(7,*) 'kin_energy_end_w:',kin_energy_end_w*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               write(7,*) 'mom_end_w:',mom_end_w*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
#endif
c
c     If KE injected onto grid is less than what it was in the dummy
c     (because of momentum cancellation due to non-static gas), add
c     the remaining energy as thermal into the host cell. add_therm_w can be greater
c     than kin_energy_dummy_w in the case where the wind slows inflowing
c     gas, because both the wind's energy and the energy of the inflow
c     are converted into heat.
c 
               add_therm_w = 0._RKIND
               kin_energy_inj_w = kin_energy_end_w - kin_energy_b_wind
               if (kin_energy_inj_w .lt. kin_energy_dummy_w) then
                  add_therm_w = kin_energy_dummy_w - kin_energy_inj_w
                  else
                     add_therm_w = 0._RKIND
                  endif
#ifdef FORTRAN_DEBUG
                  write(7,*) "canceled wind mom injected as thermal:",
     &                  add_therm_w*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
#endif
               endif
c
c           Removing the residual momentum to avoid center of mass
c           drift removes kinetic energy, so add that back in as heat
               add_therm_w = add_therm_w + e_res
c           Update energy fields to deposit cancellation + residual heat
c           into only the host cell
               te(ic,jc,kc) = te(ic,jc,kc) + 
     &                           add_therm_w/d(ic,jc,kc)
               if (idual .eq. 1_IKIND) ge(ic,jc,kc) = ge(ic,jc,kc) + 
     &                           add_therm_w/d(ic,jc,kc)
c
c           Convert momenta back to velocities and transform back to lab frame
c
               call momentum_convert(u, v, w, d, up_avg_w, vp_avg_w,
     &               wp_avg_w, nx, ny, nz, ic, jc, kc, -1_IKIND,
     &               distrad, diststep)
c
c           Add the change in the kinetic energy to the total energy field
c
               do k = -distrad, +distrad
                  stepk = abs(k)
                  do j = -distrad, +distrad
                     stepj = stepk + abs(j)
                     do i = -distrad, +distrad
                        cellstep = stepj + abs(i)
                        if (cellstep .le. diststep) then
                           vel_mag = (u(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                       v(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                       w(ic+i ,jc+j ,kc+k)**2_IKIND)**0.5_RKIND
                           ke_after_w = 0.5_RKIND *
     &                       d(ic+i, jc+j, kc+k) * vel_mag**2_IKIND
                           delta_ke_w = ke_after_w -
     &                    ke_b_wind(i+distrad+1,j+distrad+1,k+distrad+1)
                           te(ic+i,jc+j,kc+k) = te(ic+i,jc+j,kc+k) +
     &                        delta_ke_w/d(ic+i,jc+j,kc+k)
                        endif
                     enddo
                  enddo
               enddo
                  else ! end pre-SN momentum injection
c
c     If only injecting mass and metals (and not momentum) from pre-SN 
c     feedback, add them to grid
c
#ifdef FORTRAN_DEBUG
c     Sum mass, energy, and momentum before. Note mom and KE calculations
c     will be wrong here because velocities haven't been converted to
c     momentum, but we only care about mass for this debugging check
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc,
     &            mass_b_wind, kin_energy_b_wind, mom_b_wind,
     &            distrad, diststep)
               write(7,*) "m before wind:", mass_b_wind*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
c
#endif
               mass_per_cell_w = pSNmass_grid(ic,jc,kc)/distcells
               metal_per_cell_w = pSNmetal_grid(ic,jc,kc)/distcells
               call add_preSN_feedback(u, v, w, d, ge, te, metal,
     &                         nx, ny, nz, ic, jc, kc,
     &                        imetal, imulti_metals, idual, itracksrc,
     &                        metalSNII, metalSNIa, useAGB, metalAGB,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell_w, 0._RKIND,
     &                        metal_per_cell_w,
     &                        dx, x1, dunits, vunits,
     &                        distrad, diststep, e_res)
#ifdef FORTRAN_DEBUG
c     Sum mass, energy, and momentum after. Note mom and KE calculations
c     will be wrong here because velocities haven't been converted to
c     momentum, but we only care about mass for this debugging check
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc,
     &            mass_end_w, kin_energy_end_w, mom_end_w,
     &            distrad, diststep)
               write(7,*) "mass_per_cell_w:", mass_per_cell_w*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
               write(7,*) "m after wind:", mass_end_w*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
c
#endif
                  endif
               endif ! end pre-SN feedback
c
c ====================================================================
c     BEGIN SUPERNOVA FEEDBACK
c
c
c     If no supernovae in this cell, skip
               if (nsn_grid(ic,jc,kc) .eq. 0._RKIND) goto 100
c
c     Before doing any feedback, move gas in host cell into the rest
c     of the injection region. This displacement is necessary for the
c     calculation of swept-up mass to be right when determining how
c     much momentum to inject. Entrain only for discrete SN events,
c     since continuous injection will slowly redistribute gas through
c     the hydro solver and doesn't require an explicit distribution.
               f_ent = 0._RKIND
               if (use_tabfbk .eq. 1_IKIND) then
                  if (stochastic .eq. 1_IKIND) f_ent = 1._RKIND
               else
                  if (exptime_in .ge. 0._RKIND) f_ent = 1._RKIND
               endif
               call entrain_host_gas(d, u, v, w, te, ge, metal,
     &                              metalSNII, metalSNIa, metalAGB,
     &                              nx, ny, nz, ic, jc, kc, f_ent,
     &                              idual, imetal, itracksrc, useAGB,
     &                              distrad, diststep, distcells,
     &                              0_IKIND, 0._RKIND)
c
c     Compute the mass and momentum that this cell is responsible for
c     injecting
c
               if (use_tabfbk .eq. 1_IKIND) then
                  energy = nsn_ii_grid(ic,jc,kc)*ergSNII
     &                  + nsn_ia_grid(ic,jc,kc)*ergSNIa
               else
                  energy = nsn*1.e51_RKIND
               endif
               mass_per_cell = meject_grid(ic,jc,kc)/distcells
               mom_per_cell = 0._RKIND
               energy_per_cell = 0._RKIND
               E_clip = 0._RKIND
               call get_mom_per_cell(d, metal, imetal, mass_per_cell,
     &                  distcells, mom_per_cell, energy_per_cell,
     &                  nsn_grid(ic,jc,kc), energy, dx, x1, dunits,
     &                  vunits, nx, ny, nz, ic, jc, kc,
     &                  mom_mult, distrad, diststep, e_cap, E_clip)
c           mom_per_cell has units of mass*velocity/volume
#ifdef FORTRAN_DEBUG
               write(7,*) 'energy_per_cell', energy_per_cell
     &                  *dunits*(dx*x1)**3._RKIND*vunits**2._RKIND
#endif
c
c     If not injecting any momentum or energy, skip
c
               mom_inj = 0._RKIND
               energy_inj = 0._RKIND
               do k = 1_IKIND,2_IKIND*distrad+1_IKIND
                  do j = 1_IKIND,2_IKIND*distrad+1_IKIND
                     do i = 1_IKIND,2_IKIND*distrad+1_IKIND
                        mom_inj = mom_inj + mom_per_cell(i, j, k)
                        energy_inj = energy_inj +
     &                      energy_per_cell(i, j, k)
                     enddo
                  enddo
               enddo
#ifdef FORTRAN_DEBUG
               write(7,*) 'energy_inj', energy_inj
     &                  *dunits*(dx*x1)**3._RKIND*vunits**2._RKIND
#endif
               if (mom_inj .eq. 0._RKIND) then
                  if (energy_inj .eq. 0._RKIND) then
#ifdef FORTRAN_DEBUG
                     write(7,*) 'no injection for cell', ic,jc,kc
#endif
                     goto 100
                  endif
               endif
c
c     Compute the kinetic energy in the affected region
c     before momentum is added.  This
c     is needed at the end of the calculation to update the
c     total energy (te) field.
c
               ke_before = 0._RKIND
               do k = -distrad, +distrad
                  stepk = abs(k)
                  do j = -distrad, +distrad
                     stepj = stepk + abs(j)
                     do i = -distrad, +distrad
                        cellstep = stepj + abs(i)
                        if (cellstep .le. diststep) then
                           ke_before(i+distrad+1, j+distrad+1,
     &                        k+distrad+1)=
     &                         0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                                 (u(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                  v(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                                  w(ic+i ,jc+j ,kc+k)**2_IKIND)
                        endif
                     enddo
                  enddo
               enddo
#ifdef FORTRAN_DEBUG
               write(7,*) 'ke_before:',ke_before
     &                  *dunits*(dx*x1)**3._RKIND*vunits**2._RKIND
#endif
c           ke_before has units of mass*velocity^2/volume
c
c     First convert velocities to momenta and transform
c     into frame comoving with explosion-weighted average velocity of particles
c     in this cell
c
               up_avg = up_grid(ic,jc,kc)/nsn_grid(ic,jc,kc)
               vp_avg = vp_grid(ic,jc,kc)/nsn_grid(ic,jc,kc)
               wp_avg = wp_grid(ic,jc,kc)/nsn_grid(ic,jc,kc)
               call momentum_convert(u, v, w, d, up_avg, vp_avg,
     &                    wp_avg, nx, ny, nz, ic, jc, kc,
     &                    +1_IKIND, distrad, diststep)
c
c           u,v,w have units of mass*velocity/volume
c           Now u,v,w are the explosion-relative *momentum* directions
c
c     Sum mass, energy, and momentum before
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc,
     &            mass_before, kin_energy_before, mom_start,
     &            distrad, diststep)
c
#ifdef FORTRAN_DEBUG
               write(7,*) 'nsn,meject,mzeject,mzeject_ii,mzeject_ia',
     &            nsn_grid(ic,jc,kc),
     &            meject_grid(ic,jc,kc)*dunits*(dx*x1)**3_IKIND
     &              /SolarMass,
     &            mzeject_grid(ic,jc,kc)*dunits*(dx*x1)**3_IKIND
     &              /SolarMass,
     &            mzeject_ii_grid(ic,jc,kc)*dunits*(dx*x1)**3_IKIND
     &              /SolarMass,
     &            mzeject_ia_grid(ic,jc,kc)*dunits*(dx*x1)**3_IKIND
     &              /SolarMass
               write(7,*) 'up, vp, wp',
     &            up_avg*vunits/1e5_IKIND,
     &            vp_avg*vunits/1e5_IKIND,
     &            wp_avg*vunits/1e5_IKIND
               write(7,*) "m before:", mass_before*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
               write(7,*) "ke before:", kin_energy_before*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               write(7,*) 'mom_start:',mom_start*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
c           mass_before has units of mass/volume
c           kin_energy_before has units of mass*velocity^2/volume
c           mom_start has units of mass*velocity/volume
c
#endif
c     Zero dummy velocity, energy, and metallicity fields
c
               do k = -distrad, distrad
                  do j = -distrad, distrad
                     do i = -distrad, distrad
                        u1(i+distrad+1,j+distrad+1,k+distrad+1)=0._RKIND
                        v1(i+distrad+1,j+distrad+1,k+distrad+1)=0._RKIND
                        w1(i+distrad+1,j+distrad+1,k+distrad+1)=0._RKIND
                        d1(i+distrad+1,j+distrad+1,k+distrad+1) =
     &                     d(i+ic,j+jc,k+kc)
                        ge1(i+distrad+1,j+distrad+1,k+distrad+1) = 
     &                     0._RKIND
                        te1(i+distrad+1,j+distrad+1,k+distrad+1) = 
     &                     0._RKIND
                        metal1(i+distrad+1,j+distrad+1,k+distrad+1) = 
     &                     0._RKIND
                     enddo
                  enddo
               enddo
c
c     Now add mass and momentum terms to local dummy fields
c     Metals and thermal energy zeroed out because we don't care about them in the dummy
c
               call add_feedback_mech(u1, v1, w1, d1, ge1, te1, metal1,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        distrad+1_IKIND, distrad+1_IKIND,
     &                        distrad+1_IKIND,
     &                        0_IKIND, imulti_metals, idual,
     &                        itracksrc,
     &                        0._RKIND, 0._RKIND, 0._RKIND,
     &                        metalSNII, metalSNIa, useAGB, metalAGB,
     &                        mass_per_cell, mom_per_cell,
     &                        dx, x1, dunits, vunits,
     &                        distrad, diststep, e_res)
c
c     Calculate mass and energy in dummy fields
c
               call sum_energy_mom(u1, v1, w1, d1,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        2_IKIND*distrad+1_IKIND,
     &                        distrad+1_IKIND, distrad+1_IKIND,
     &                        distrad+1_IKIND,
     &                        mass_dummy, kin_energy_dummy, mom_dummy,
     &                        distrad, diststep)
#ifdef FORTRAN_DEBUG
               write(7,*) "m dummy:", mass_dummy*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
               write(7,*) "ke dummy:", kin_energy_dummy*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               write(7,*) 'mom_dummy:',mom_dummy*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
#endif
c
c     Now deposit the momentum, mass, and metals onto the actual grid.
c
               metal_per_cell = mzeject_grid(ic,jc,kc) / distcells
               metal_ii_per_cell = mzeject_ii_grid(ic,jc,kc)/distcells
               metal_ia_per_cell = mzeject_ia_grid(ic,jc,kc)/distcells
               call add_feedback_mech(u, v, w, d, ge, te, metal,
     &                         nx, ny, nz, ic, jc, kc,
     &                        imetal, imulti_metals, idual,
     &                        itracksrc,
     &                        metal_per_cell,
     &                        metal_ia_per_cell,
     &                        metal_ii_per_cell,
     &                        metalSNII, metalSNIa, useAGB, metalAGB,
     &                        mass_per_cell, mom_per_cell,
     &                        dx, x1, dunits, vunits,
     &                        distrad, diststep, e_res)
c
c     Calculate the mass and energy in the affected region again.
c
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &                          ic, jc, kc,
     &                          mass_end, kin_energy_end, mom_end,
     &                          distrad, diststep)
#ifdef FORTRAN_DEBUG
               write(7,*) "m end:", mass_end*
     &               (dunits*(dx*x1)**3._RKIND)/SolarMass
               write(7,*) 'kin_energy_end:',kin_energy_end*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               write(7,*) 'mom_end:',mom_end*dunits*vunits
     &                  *(dx*x1)**3._RKIND/SolarMass/1e5_RKIND
#endif
c
c     If KE injected onto grid is less than what it was in the dummy
c     (because of momentum cancellation due to non-static gas), calculate
c     how much energy to inject as thermal. add_therm can be greater
c     than kin_energy_dummy in the case where an explosion slows inflowing
c     gas, because both the explosion's energy and the energy of the inflow
c     are converted into heat.
c
               kin_energy_inj = kin_energy_end - kin_energy_before
               if (kin_energy_inj .lt. kin_energy_dummy) then
                  add_therm = kin_energy_dummy - kin_energy_inj
#ifdef FORTRAN_DEBUG
                  write(7,*) "leftover energy injected as thermal:",
     &                  add_therm*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
#endif
c
               else
                  add_therm = 0._RKIND
               endif ! end lost kinetic -> thermal energy injection
c
c           Removing the residual momentum to avoid center of mass
c           drift removes kinetic energy, so add that back in as heat
               add_therm = add_therm + e_res
c
c           Convert momenta back to velocities and transform back to lab frame
c
               call momentum_convert(u, v, w, d, up_avg, vp_avg,
     &               wp_avg, nx, ny, nz, ic, jc, kc, -1_IKIND,
     &               distrad, diststep)
c
c           Add the change in the kinetic energy to the total energy field,
c           and inject energy lost due to colliding flows (add_therm) and
c           thermal energy from resolved Sedov-Taylor phase (energy_per_cell).
c
               ke_injected = 0._RKIND
               te_sum = 0._RKIND
               ge_sum = 0._RKIND
               ke_sum = 0._RKIND
               mom_sum = 0._RKIND
               mass_sum = 0._RKIND
               delta_ge = add_therm/mass_end
               do k = -distrad, +distrad
                  stepk = abs(k)
                  do j = -distrad, +distrad
                     stepj = stepk + abs(j)
                     do i = -distrad, +distrad
                        cellstep = stepj + abs(i)
                        if (cellstep .le. diststep) then
                        add_ge = delta_ge +
     &                     energy_per_cell(i+distrad+1,
     &                     j+distrad+1,k+distrad+1)/d(ic+i,jc+j,kc+k)
                        te(ic+i,jc+j,kc+k) = te(ic+i,jc+j,kc+k) + add_ge
                        if (idual .eq. 1_IKIND) then
                           ge(ic+i,jc+j,kc+k) = ge(ic+i,jc+j,kc+k) +
     &                       add_ge
                           ge_sum = ge_sum +
     &                       ge(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k)
                        endif
                        vel_mag = (u(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                    v(ic+i ,jc+j ,kc+k)**2_IKIND +
     &                    w(ic+i ,jc+j ,kc+k)**2_IKIND)**0.5_RKIND
                        ke_after = 0.5_RKIND * d(ic+i, jc+j, kc+k) *
     &                    vel_mag**2_IKIND
                        delta_ke = ke_after -
     &                    ke_before(i+distrad+1,j+distrad+1,k+distrad+1)
                        te(ic+i,jc+j,kc+k) = te(ic+i,jc+j,kc+k) +
     &                     delta_ke/d(ic+i,jc+j,kc+k)
                        ke_injected = ke_injected + delta_ke
                        mass_term = d(ic+i,jc+j,kc+k)
                        ke_sum = ke_sum + ke_after
                        mom_sum = mom_sum + mass_term*vel_mag
                        mass_sum = mass_sum + mass_term
                        te_sum = te_sum +
     &                    te(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k)
c                       ke_after has units of mass*velocity^2/volume
c                       te has units of mass*velocity^2/volume/mass*volume = velocity^2
c                       te is *specific* energy (energy/mass)
                        endif
                     enddo
                  enddo
               enddo
#ifdef FORTRAN_DEBUG
               write(7,*) 'mass, mom, KE_sum, KE_inj, GE, TE:',
     &          mass_sum*dunits*(x1*dx)**3._RKIND/SolarMass,
     &          mom_sum*dunits*(x1*dx)**3._RKIND*vunits/
     &          SolarMass/1e5_RKIND,
     &          ke_sum*dunits*(x1*dx)**3._RKIND*vunits**2._RKIND,
     &          ke_injected*dunits*(x1*dx)**3._RKIND*vunits**2._RKIND,
     &          ge_sum*dunits*(x1*dx)**3._RKIND*vunits**2._RKIND,
     &          te_sum*dunits*(x1*dx)**3._RKIND*vunits**2._RKIND
#endif
c
               if (feedback_log .eq. 1_IKIND) then
                  write(5,*) t*t1/yr_s/1e6_RKIND, ic, jc, kc,
     &             nsn_grid(ic,jc,kc), meject_grid(ic,jc,kc)*
     &             (dunits*(x1*dx)**3._RKIND)/SolarMass,
     &             mzeject_grid(ic,jc,kc)*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass,
     &             mom_inj*(dunits*(dx*x1)**3._RKIND)*
     &             vunits/SolarMass/1e5_RKIND,
     &             energy_inj*dunits*
     &             vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND,
     &             add_therm*dunits*
     &             vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND,
     &             kin_energy_dummy*dunits*
     &             vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND,
     &             kin_energy_inj*dunits*
     &             vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND,
     &             pSNmass_grid(ic,jc,kc)*
     &             (dunits*(x1*dx)**3._RKIND)/SolarMass,
     &             pSNmetal_grid(ic,jc,kc)*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass,
     &             pSNmom_grid(ic,jc,kc)*(dunits*(dx*x1)**3._RKIND)*
     &             vunits/SolarMass/1e5_RKIND,
     &             add_therm_w*dunits*
     &             vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               endif
c
 100           continue
c
c ====================================================================
c     BEGIN AGB WIND FEEDBACK
c
c     Deposit AGB wind mass and metals directly into the star
c     particles' host cell, with no momentum or energy injection.
c     This is done after SNe (and pre-SN) feedback so that it does
c     not interfere with the mass/energy/momentum bookkeeping above.
c
               if (useAGB .eq. 1_IKIND .and.
     &             agbmass_grid(ic,jc,kc) .gt. 0._RKIND) then
                  if (imetal .eq. 1_IKIND) then
                     metal(ic,jc,kc) = (metal(ic,jc,kc)*d(ic,jc,kc) +
     &                 agbmetal_grid(ic,jc,kc))
     &                     / (d(ic,jc,kc) + agbmass_grid(ic,jc,kc))
                     if (metal(ic,jc,kc) .gt. 0.95_RKIND)
     &                  metal(ic,jc,kc) = 0.95_RKIND
                     if (itracksrc .eq. 1_IKIND) then
                        metalAGB(ic,jc,kc) =
     &                    (metalAGB(ic,jc,kc)*d(ic,jc,kc) +
     &                    agbmetal_grid(ic,jc,kc))
     &                    / (d(ic,jc,kc) + agbmass_grid(ic,jc,kc))
c              Dilute other metal fields by the amount of mass added
c              by AGB winds so they don't gain metals from AGB
                        metalSNII(ic,jc,kc) =
     &                        (metalSNII(ic,jc,kc)*d(ic,jc,kc))
     &                        / (d(ic,jc,kc) + agbmass_grid(ic,jc,kc))
                        metalSNIa(ic,jc,kc) =
     &                        (metalSNIa(ic,jc,kc)*d(ic,jc,kc))
     &                        / (d(ic,jc,kc) + agbmass_grid(ic,jc,kc))
                        if (metalAGB(ic,jc,kc) .gt. 0.95_RKIND)
     &                    metalAGB(ic,jc,kc) = 0.95_RKIND
                     endif
                     endif
                     d(ic,jc,kc) = d(ic,jc,kc) + agbmass_grid(ic,jc,kc)
               endif ! end AGB wind feedback
c
            enddo
         enddo
      enddo ! end loop over feedback grid
c
c
      if (feedback_log .eq. 1) then
         close(4)
         close(5)
      endif
#ifdef FORTRAN_DEBUG
      close(7)
      close(8)
#endif
      return
      end
c
c
c ==========================================================
c
c     Entrain f_ent*(1-beta) of the host cell's gas into the other
c     cells of the injection region (Kimm & Cen 2014, eq. A6), with
c     beta = 1/distcells, so each receiving cell gets f_ent*d_c/
c     distcells. Conserves mass, momentum, total energy, metals and
c     every passive field. Works on lab-frame velocities, so it must
c     be called before momentum_convert(+1). The host cell keeps its
c     specific quantities (velocity, te, ge, fractions); only its
c     density (and any density-form passive fields) drop.
c
      subroutine entrain_host_gas(d, u, v, w, te, ge, metal,
     &     metalSNII, metalSNIa, metalAGB,
     &     nx, ny, nz, ic, jc, kc, f_ent, idual, imetal, itracksrc,
     &     useAGB, distrad, diststep, distcells, distmode, sphererad2)
c
      implicit none
#include "fortran_types.def"
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual, imetal, itracksrc,
     &          useAGB, distrad, diststep, distcells, distmode
      R_PREC    d(nx,ny,nz), u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),
     &          te(nx,ny,nz), ge(nx,ny,nz), metal(nx,ny,nz),
     &          metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz),
     &          metalAGB(nx,ny,nz)
      R_PREC    f_ent, sphererad2
c
      INTG_PREC i, j, k, n, ii, jj, kk
      LOGIC_PREC inregion
      R_PREC    dc, uc, vc, wc, tec, gec, Zc, ZIIc, ZIac, ZAGBc,
     &          dmove, xmove, xkeep, dold, dnew, ke_old, ke_new
c
      if (f_ent .le. 0._RKIND .or. distcells .le. 1_IKIND) return
c
      dc  = d(ic,jc,kc)
      uc  = u(ic,jc,kc)
      vc  = v(ic,jc,kc)
      wc  = w(ic,jc,kc)
      tec = te(ic,jc,kc)
      gec = ge(ic,jc,kc)
      Zc  = 0._RKIND
      ZIIc = 0._RKIND
      ZIac = 0._RKIND
      ZAGBc = 0._RKIND
      if (imetal .eq. 1_IKIND) Zc = metal(ic,jc,kc)
      if (imetal .eq. 1_IKIND .and. itracksrc .eq. 1_IKIND) then
         ZIIc = metalSNII(ic,jc,kc)
         ZIac = metalSNIa(ic,jc,kc)
         if (useAGB .eq. 1_IKIND) ZAGBc = metalAGB(ic,jc,kc)
      endif
c
c     Density moved into each receiving cell, and the fraction of the
c     host's gas that this is
      dmove = f_ent*dc/real(distcells, RKIND)
      xmove = dmove/dc
      xkeep = 1._RKIND - real(distcells-1, RKIND)*xmove
c
      do k = -distrad, distrad
         do j = -distrad, distrad
            do i = -distrad, distrad
               if (distmode .eq. 1_IKIND) then
                  inregion = real(i*i+j*j+k*k, RKIND) .le. sphererad2
               else
                  inregion = abs(i)+abs(j)+abs(k) .le. diststep
               endif
c     If in injection region AND not the center cell:
               if (inregion .and.
     &           (i .ne. 0 .or. j .ne. 0 .or. k .ne. 0)) then
                  dold = d(ic+i,jc+j,kc+k)
                  dnew = dold + dmove
c     This cell's KE is its original KE + KE of gas displaced into this
c     cell, which carries the velocity of the host cell it came from
                  ke_old = 0.5_RKIND*dold*(u(ic+i,jc+j,kc+k)**2
     &                     + v(ic+i,jc+j,kc+k)**2
     &                     + w(ic+i,jc+j,kc+k)**2)
     &                     + 0.5_RKIND*dmove*(uc**2 + vc**2 + wc**2)
c     Update the velocity by assuming conservation of momentum from
c     moving material from one cell to another. Note this is not a
c     momentum kick, just purely displacement of mass that had a velocity.
                  u(ic+i,jc+j,kc+k) = (dold*u(ic+i,jc+j,kc+k)
     &                                + dmove*uc)/dnew
                  v(ic+i,jc+j,kc+k) = (dold*v(ic+i,jc+j,kc+k)
     &                                + dmove*vc)/dnew
                  w(ic+i,jc+j,kc+k) = (dold*w(ic+i,jc+j,kc+k)
     &                                + dmove*wc)/dnew
c     This cell's new KE is calculated from its new density and velocity,
c     not guaranteed to be equal to the sum of the KEs of the original
c     cell and the moving material
                  ke_new = 0.5_RKIND*dnew*(u(ic+i,jc+j,kc+k)**2
     &                     + v(ic+i,jc+j,kc+k)**2
     &                     + w(ic+i,jc+j,kc+k)**2)
c     Update total energy by the amount of total energy that has been
c     moved into this cell
                  te(ic+i,jc+j,kc+k) = (dold*te(ic+i,jc+j,kc+k)
     &                                 + dmove*tec)/dnew
                  if (idual .eq. 1_IKIND)
c     If the new KE is different by more than just the fraction of KE
c     that has been moved into this cell, add the difference as heat
     &               ge(ic+i,jc+j,kc+k) = (dold*ge(ic+i,jc+j,kc+k)
     &                                    + dmove*gec
     &                                    + (ke_old - ke_new))/dnew
c     Update all metal fractions as a mass-weighted mix
                  if (imetal .eq. 1_IKIND) then
                     metal(ic+i,jc+j,kc+k) = (dold*metal(ic+i,jc+j,kc+k)
     &                                       + dmove*Zc)/dnew
                     if (itracksrc .eq. 1_IKIND) then
                        metalSNII(ic+i,jc+j,kc+k) =
     &                            (dold*metalSNII(ic+i,jc+j,kc+k)
     &                            + dmove*ZIIc)/dnew
                        metalSNIa(ic+i,jc+j,kc+k) =
     &                            (dold*metalSNIa(ic+i,jc+j,kc+k)
     &                            + dmove*ZIac)/dnew
                        if (useAGB .eq. 1_IKIND)
     &                     metalAGB(ic+i,jc+j,kc+k) =
     &                              (dold*metalAGB(ic+i,jc+j,kc+k)
     &                              + dmove*ZAGBc)/dnew
                     endif
                  endif
c     Finally update the density in this cell itself
                  d(ic+i,jc+j,kc+k) = dnew
               endif ! end in-region
            enddo
         enddo
      enddo
c
c     Update host cell with the kept fraction of density.
c     Energies and metallicities don't need updating because they're
c     specific energies and metal fractions so already scale with density
      d(ic,jc,kc) = xkeep*dc
c
      return
      end
c
c ==========================================================
c
c     Calculate the amount of momentum to add to each neighboring cell
c     Formulae taken from Kimm & Cen (2014), Appendix A
c
      subroutine get_mom_per_cell(d, metal, imetal, mass_per_cell,
     &               distcells, mom_per_cell, energy_per_cell,
     &               nsn, E_tot, dx, x1, dunits, vunits,
     &               nx, ny, nz, ic, jc, kc,
     &               mom_mult, distrad, diststep, e_cap, E_clip)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, distrad, diststep, distcells
      INTG_PREC imetal
      R_PREC    d(nx,ny,nz), metal(nx,ny,nz), dx
      R_PREC    mom_per_cell(2*distrad+1, 2*distrad+1, 2*distrad+1)
      R_PREC    energy_per_cell(2*distrad+1, 2*distrad+1, 2*distrad+1)
      R_PREC    mass_per_cell, nsn, E_tot, x1, dunits, vunits
      R_PREC    mom_mult, e_cap, E_clip
c
c     Locals
c
      INTG_PREC i, j, k, dist, cellstep, stepk, stepj
      LOGIC_PREC inregion
      R_PREC    Nnbors, cell_vol, ecode, Mej, dMej, E_nb, frac_sn,
     &          E_ev51, den_cell, Z_cell, Zp, n_H, pSN_term, chi_th,
     &          dMswept, chi, fe, pSN, ke_from_p, s_cap, E_cell
      R_PREC    Xh

      Xh = 0.76_RKIND
c
      Nnbors = real(distcells, RKIND) - 1._RKIND ! 26
      cell_vol = (dx*x1)**3
      ecode = dunits*cell_vol*vunits**2 ! cell energy density units
      Mej = mass_per_cell*real(distcells, RKIND) * dunits*cell_vol
      dMej = Mej/real(distcells, RKIND)
c
      E_nb = E_tot/Nnbors
c     N >= 1 is one event of N SNe; N < 1 is a fraction of one
      frac_sn = min(nsn, 1._RKIND)
      E_ev51 = E_tot/frac_sn/1.e51_RKIND
c
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
               if (cellstep .le. diststep .and. dist .gt. 0._IKIND) then
c
c        If not using momentum at all through mom_mult = 0, set energy
c        per cell as simply split between all cells in injection zone
                  if (mom_mult .eq. 0._RKIND) then
                     energy_per_cell(i+distrad+1,
     &                    j+distrad+1,k+distrad+1) = E_nb/ecode
                  else
c
c       Compute the mass fraction of swept-up to ejected mass.
c
                  den_cell = d(ic+i,jc+j,kc+k)*dunits
                  Z_cell = 0._RKIND
                  if (imetal .eq. 1_IKIND)
     &                        Z_cell = metal(ic+i,jc+j,kc+k)
                  if (Z_cell .ne. Z_cell) Z_cell = 0._RKIND
                  Z_cell = max(0._RKIND, min(Z_cell, 0.95_RKIND))
c        Scheme requires a floor on metallicity
                  Zp  = max(Z_cell/0.02_RKIND, 0.01_RKIND)
c        Derive H density from metallicity, assuming Solar relative abundances
                  n_H = max(Xh*(1._RKIND - Z_cell)*den_cell/mass_h,
     &                   1.e-6_RKIND)
c        Define terminal momentum per cell and threshold for transition
c        between energy-conserving and momentum-conserving phases
c        Kimm & Cen (2014) eq. A4
                  pSN_term = frac_sn*3.e10_RKIND*SolarMass*
     &                     E_ev51**(16._RKIND/17._RKIND)*
     &                     n_H**(-2._RKIND/17._RKIND)*Zp**(-0.14_RKIND)
c        The energy-conserving and momentum-conserving branches meet
c        at chi_th, calculated by setting pSN_term = pSN defined below.
c        When chi = chi_th, fe = 2/3, and the additional factor of 2
c        in pSN below means solving for chi_th gives a factor of 0.75
                  chi_th = 0.75_RKIND*pSN_term**2/(Mej*E_tot)
c        chi_th can get very small if a small fraction of SNe are going
c        off (in continuous case) so floor it just above 1
                  chi_th = max(chi_th, 1._RKIND + 1.e-6_RKIND)
c
c        Swept-up mass is the sum of the mass this cell already had,
c        which includes already entrained mass displaced from the
c        supernova host cell, and its share of the ejecta from the SN
c
                  dMswept = den_cell*cell_vol + dMej
                  chi = dMswept/dMej
c
c        Compare chi to threshold and set momentum injection.
c
                  if (chi .lt. chi_th) then ! energy-conserving, resolved
c        Momentum comes from eq. A1
                     fe = 1._RKIND - ((chi - 1._RKIND)/
     &                        (3._RKIND*(chi_th - 1._RKIND)))
                     pSN = mom_mult * sqrt(2._RKIND*chi*Mej*fe*E_tot)
c
c        This is the resolved case, so need to inject thermal energy
c        to get up to E_tot per SNe. Use pSN to calculate kinetic energy
c        and inject thermal into this cell equal to total - kinetic
                     ke_from_p = 0.5_RKIND*(pSN/Nnbors)**2./dMswept
                     energy_per_cell(i+distrad+1,
     &                               j+distrad+1,
     &                               k+distrad+1) =
     &                  max(E_nb - ke_from_p, 0._RKIND)/ecode
                  else
c        Unresolved, only momentum
                     pSN = mom_mult * pSN_term
                  endif                     
c                 
c        Cap the specific energy deposited in this neighbor (kinetic +
c        thermal) at e_cap [erg/g], which bounds both the kick speed and
c        the post-SN sound speed. Momentum scales by sqrt(s) and heat by
c        s, so the kinetic/thermal split of this cell is unchanged.
c        Clipped energy is discarded and tallied in E_clip [erg].
                  if (e_cap .gt. 0._RKIND) then
                     E_cell = 0.5_RKIND*(pSN/Nnbors)**2/dMswept +
     &                  energy_per_cell(i+distrad+1,j+distrad+1,
     &                                  k+distrad+1)*ecode
                     if (E_cell .gt. e_cap*dMswept) then
                        s_cap = e_cap*dMswept/E_cell
                        pSN = pSN*sqrt(s_cap)
                        energy_per_cell(i+distrad+1,j+distrad+1,
     &                     k+distrad+1) = s_cap*
     &                     energy_per_cell(i+distrad+1,j+distrad+1,
     &                                     k+distrad+1)
                        E_clip = E_clip + (1._RKIND - s_cap)*E_cell
                     endif
                  endif
c        pSN has units of cgs-mass*cgs-velocity, want it to have units of code-mass*code-velocity/code-volume
c        conversion of mass from code units is (dunits*(dx*x1)**3._RKIND) (density times volume)
                  pSN = pSN/(dunits*cell_vol)/vunits/Nnbors
c        dunits = cgs-density/code-density = cgs-mass*code-volume/(cgs-volume*code-mass)
c        x1^3 = cgs-volume/code-volume
c        vunits = cgs-velocity/code-velocity
c        dx^3 = code-volume
c        dunits*(dx*x1)**3._RKIND = cgs-mass*code-volume/(cgs-volume*code-mass)*code-volume*cgs-volume/code-volume
c        dunits*(dx*x1)**3._RKIND = cgs-mass*code-volume/code-mass
c        pSN = cgs-mass*cgs-velocity/(cgs-mass*code-volume/code-mass)/(cgs-velocity/code-velocity)
c        pSN = code-mass*code-velocity/code-volume
c        pSN now has units code-mass*code-velocity/code-volume
c        So pSN is actually in the right units after all!
c
                  mom_per_cell(i+distrad+1,j+distrad+1,k+distrad+1)=pSN
c
                  endif
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Convert velocities to momentum and back
c
      subroutine momentum_convert(u, v, w, d, up, vp, wp,
     &                    nx, ny, nz, ic, jc, kc, idir, 
     &                    distrad, diststep)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, distrad, diststep
      INTG_PREC idir
      R_PREC    d(nx, ny, nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    up, vp, wp
c
c     Locals
c
      INTG_PREC i, j, k, cellstep, stepk, stepj
c
c     Error check
c
      if (idir .ne. -1_IKIND .and. idir .ne. 1_IKIND) then
         write(6,*) 'incorrect idir value in momentum call', idir
         ERROR_MESSAGE
      endif
c
c     Loop over velocities, multiplying by densities (or dividing if
c       converting back)
c
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               if (cellstep .le. diststep) then
c
c              idir = +1: convert vel -> mom
c
               if (idir. eq. +1_IKIND) then
                  u(ic+i ,jc+j ,kc+k) = (u(ic+i ,jc+j ,kc+k)-up) *
     &                       d(ic+i, jc+j, kc+k)
                  v(ic+i ,jc+j ,kc+k) = (v(ic+i ,jc+j ,kc+k)-vp) *
     &                       d(ic+i, jc+j, kc+k)
                  w(ic+i ,jc+j ,kc+k) = (w(ic+i ,jc+j ,kc+k)-wp) *
     &                       d(ic+i, jc+j, kc+k)
c
c              if idir = -1: convert mom -> vel
c
               else
                  u(ic+i ,jc+j ,kc+k) = u(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + up
                  v(ic+i ,jc+j ,kc+k) = v(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + vp
                  w(ic+i ,jc+j ,kc+k) = w(ic+i ,jc+j ,kc+k) /
     &                    (d(ic+i, jc+j, kc+k)) + wp
               endif
c
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy
c       Note that pu, pv, pw are momenta
c
      subroutine sum_energy_mom(pu, pv, pw, d,
     &                           nx, ny, nz,
     &                           ic, jc, kc,
     &                           mass_sum, kin_energy_sum, mom_sum,
     &                           distrad, diststep)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, distrad, diststep
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz)
      R_PREC    mass_sum, kin_energy_sum
      R_PREC    mom_sum
c
c     Locals
c
      INTG_PREC i, j, k, cellstep, stepk, stepj
      R_PREC    mass_term, mom_term, kin_energy
c
      mass_sum = 0._RKIND
      kin_energy_sum = 0._RKIND
      mom_sum = 0._RKIND
c
c     Sum mass and energy      
c
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               if (cellstep .le. diststep) then
               mass_term = d(ic+i, jc+j, kc+k)
               mom_term = pu(ic+i, jc+j, kc+k)**2_IKIND +
     &                    pv(ic+i, jc+j, kc+k)**2_IKIND +
     &                    pw(ic+i, jc+j, kc+k)**2_IKIND
c              mass_term has units of mass/volume
c              mom_term has units of mass^2*velocity^2/volume^2
c
c              Compute total mass and energy
c              
               kin_energy = mom_term / (2._RKIND * mass_term)
c              kin_energy has units of mass^2*velocity^2/volume^2 / (mass/volume)
c                          = mass*velocity^2/volume
               mass_sum = mass_sum + mass_term
               kin_energy_sum = kin_energy_sum + kin_energy
               mom_sum = mom_sum + mom_term**0.5_RKIND
c              mom_grid has units of mass*velocity/volume
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Add mass, momentum, and energy to affected cells.
c       Note that pu, pv, pw are momenta
c
      subroutine add_feedback_mech(pu, pv, pw, d, ge, te, metal,
     &                        nx, ny, nz, ic, jc, kc,
     &                        imetal, imulti_metals, idual,
     &                        itracksrc,
     &                        mzeject, mzeject_ia, mzeject_ii,
     &                        metalSNII, metalSNIa, useAGB, metalAGB,
     &                        mass_per_cell, mom_per_cell,
     &                        dx, x1, dunits, vunits,
     &                        distrad, diststep, e_res)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual, dist, distrad, diststep
      INTG_PREC imetal, imulti_metals, useAGB
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_per_cell, e_res
      R_PREC    mom_per_cell(2*distrad+1,2*distrad+1,2*distrad+1)
      R_PREC    mzeject, mzeject_ia, mzeject_ii
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    metalAGB(nx,ny,nz)
      R_PREC    dunits, x1, vunits, dx
      R_PREC    p_start, p_end
      INTG_PREC    itracksrc
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1, cellstep, stepk, stepj
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw,
     &          delta_p, vel_kick, dratio, r
      R_PREC    pu_start, pu_deposit, pv_start, pv_deposit, pw_start,
     &          pw_deposit
      R_PREC    p_res(3), m_res, f_res, pk
c
c     Calculate the net momentum all the radial kicks will deposit.
c     This will need to be subtracted off of every cell in the feedback
c     zone to ensure the net momentum of the region stays zero (we
c     don't want the whole region to drift).
c     The host cell gets no radial kick but its mass is part of the system
c     whose center of mass must stay at rest in the star frame.
c
      p_res(1) = 0._RKIND
      p_res(2) = 0._RKIND
      p_res(3) = 0._RKIND
      m_res = 0._RKIND
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               if (cellstep .le. diststep) then
                  m_res = m_res + d(ic+i,jc+j,kc+k) + mass_per_cell
                  dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
                  if (dist .gt. 0_IKIND) then
                     pk = mom_per_cell(i+distrad+1,j+distrad+1,
     &                                 k+distrad+1)/
     &                    sqrt(real(dist, RKIND))
                     p_res(1) = p_res(1) + real(i, RKIND)*pk
                     p_res(2) = p_res(2) + real(j, RKIND)*pk
                     p_res(3) = p_res(3) + real(k, RKIND)*pk
                  endif
               endif
            enddo
         enddo
      enddo
c     Removing a bulk momentum drift would remove some kinetic energy
c     that should be deposited, so sum up how much and deposit it
c     as thermal energy later
      e_res = 0.5_RKIND*(p_res(1)**2 + p_res(2)**2 + p_res(3)**2)/m_res
c
c     Loop cells in particle-frame
c
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               if (cellstep .le. diststep) then
c              Calculate distance of this cell from host cell
               dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
c              Decompose the radial momentum onto Cartesian axes using
c              the outward unit vector (i,j,k)/r where r = ||(i,j,k)||.
c              At the center cell i=j=k=0 so delta_p components are zero;
c              r=1 avoids division by zero in that case.
               if (dist .eq. 0_IKIND) then
                  r = 1._RKIND
               else
                  r = sqrt(real(dist, RKIND))
               endif
               delta_mass = mass_per_cell
c              delta_mass has units of mass/volume
c              Remove this cell's mass-weighted share of the bulk momentum
               f_res = (d(ic+i,jc+j,kc+k) + delta_mass)/m_res
               delta_pu = real(i, RKIND)/r *
     &            mom_per_cell(i+distrad+1,j+distrad+1,k+distrad+1)
     &            - f_res*p_res(1)
               delta_pv = real(j, RKIND)/r *
     &            mom_per_cell(i+distrad+1,j+distrad+1,k+distrad+1)
     &            - f_res*p_res(2)
               delta_pw = real(k, RKIND)/r *
     &            mom_per_cell(i+distrad+1,j+distrad+1,k+distrad+1)
     &            - f_res*p_res(3)
               delta_p = (delta_pu**2._RKIND + delta_pv**2._RKIND + 
     &            delta_pw**2._RKIND)**(0.5_RKIND)
#ifdef FORTRAN_DEBUG
               write(7,*) 'i,j,k, dist', i, j, k, dist
               write(7,*) 'cellstep, diststep', cellstep, diststep
               write(7,*) 'p_before:',pu(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            pv(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            pw(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND
               write(7,*) 'delta_p:',delta_pu*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            delta_pv*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            delta_pw*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND
#endif
c
c     Add momentum
c
               pu_start = pu(ic+i,jc+j,kc+k)
               pv_start = pv(ic+i,jc+j,kc+k)
               pw_start = pw(ic+i,jc+j,kc+k)
               p_start = (pu_start**2._RKIND + pv_start**2._RKIND +
     &           pw_start**2.)**(0.5_RKIND)
               pu(ic+i,jc+j,kc+k) =
     &                       pu_start + delta_pu
               pv(ic+i,jc+j,kc+k) =
     &                       pv_start + delta_pv
               pw(ic+i,jc+j,kc+k) =
     &                       pw_start + delta_pw
               p_end = (pu(ic+i,jc+j,kc+k)**2._RKIND +
     &                  pv(ic+i,jc+j,kc+k)**2._RKIND +
     &                  pw(ic+i,jc+j,kc+k)**2._RKIND)**(0.5_RKIND)
#ifdef FORTRAN_DEBUG
               write(7,*) 'v_start, v_final, mass:',
     &            p_start/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND,
     &            p_end/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND,
     &            (d(ic+i,jc+j,kc+k) + delta_mass)
     &              *dunits*(dx*x1)**3_IKIND/SolarMass
#endif
c
c     Adjust thermal energy to account for new mass added
c
               dratio = d(i+ic,j+jc,k+kc)/
     &           (d(ic+i,jc+j,kc+k) + delta_mass)
               te(ic+i,jc+j,kc+k) = te(i+ic,j+jc,k+kc)*dratio
               if (idual .eq. 1_IKIND)
     &           ge(ic+i,jc+j,kc+k) = ge(ic+i,jc+j,kc+k)*dratio
c
c     Metal feedback
c     metal fields are metal fractions here, need to convert to metal density
c     before adding new metals
c
               if (imetal .eq. 1_IKIND) then
                  metal(ic+i,jc+j,kc+k) =
     &               (metal(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) +mzeject)
     &               /(d(ic+i,jc+j,kc+k) + delta_mass)
c                 Ensure the metallicity doesn't get too high
                  if (metal(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &              metal(ic+i,jc+j,kc+k) = 0.95_RKIND
                  if (itracksrc .eq. 1) then
                     metalSNII(ic+i,jc+j,kc+k) =
     &                 (metalSNII(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) +
     &                 mzeject_ii)/(d(ic+i,jc+j,kc+k) + delta_mass)
c                 Put same cap on other metal fields
                     if (metalSNII(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &                 metalSNII(ic+i,jc+j,kc+k) = 0.95_RKIND
                     metalSNIa(ic+i,jc+j,kc+k) =
     &                 (metalSNIa(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) +
     &                 mzeject_ia)/(d(ic+i,jc+j,kc+k) + delta_mass)
                     if (metalSNIa(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &                 metalSNIa(ic+i,jc+j,kc+k) = 0.95_RKIND
c                 Dilute AGB metal field by the amount of non-AGB mass
                     if (useAGB .eq. 1_IKIND) metalAGB(ic+i,jc+j,kc+k) = 
     &                 (metalAGB(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k))
     &                 / (d(ic+i,jc+j,kc+k) + delta_mass)
                  endif
               endif
c
c     Add mass
c
               d(ic+i,jc+j,kc+k) = d(ic+i,jc+j,kc+k)
     &                       + delta_mass
c
c     End loop over cells in particle-frame
c
               endif
            enddo
         enddo
      enddo
      return
      end
c
c ==========================================================
c
c     Add mass, metals, and momentum from pre-SN feedback
c       Note that pu, pv, pw are momenta
c
      subroutine add_preSN_feedback(pu, pv, pw, d, ge, te, metal,
     &                        nx, ny, nz, ic, jc, kc,
     &                        imetal, imulti_metals, idual, itracksrc,
     &                        metalSNII, metalSNIa, useAGB, metalAGB,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell, mom_per_cell,
     &                        metal_per_cell,
     &                        dx, x1, dunits, vunits,
     &                        distrad, diststep, e_res)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual, dist, distrad, diststep
      INTG_PREC imetal, imulti_metals, kick_cap, preSNmom, useAGB,
      INTG_PREC itracksrc
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    metalAGB(nx,ny,nz)
      R_PREC    mass_per_cell, mom_per_cell, e_res
      R_PREC    metal_per_cell
      R_PREC    dunits, x1, vunits, dx
      R_PREC    p_start, p_end
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1, cellstep, stepk, stepj
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw,
     &          delta_p, vel_kick, dratio, r
      R_PREC    pu_start, pu_deposit, pv_start, pv_deposit, pw_start,
     &          pw_deposit
      R_PREC    p_res(3), m_res, f_res, pk, m_nb, f_nb
c
c     First calculate total mass in region
      m_res = 0._RKIND
      m_nb = 0._RKIND
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               if (cellstep .le. diststep) then
                  dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
                  m_res = m_res + d(ic+i,jc+j,kc+k) + mass_per_cell
                  if (dist .gt. 0_IKIND)
     &               m_nb = m_nb + d(ic+i,jc+j,kc+k) + mass_per_cell
               endif
            enddo
         enddo
      enddo
c
c     Calculate the net momentum all the radial kicks will deposit.
c     This will need to be subtracted off of every cell in the feedback
c     zone to ensure the net momentum of the region stays zero (we
c     don't want the whole region to drift).
c     The host cell gets no radial kick but its mass is part of the system
c     whose center of mass must stay at rest in the star frame.
c
      p_res(1) = 0._RKIND
      p_res(2) = 0._RKIND
      p_res(3) = 0._RKIND
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               if (cellstep .le. diststep) then
                  dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
                  if (dist .gt. 0_IKIND) then
c              normalize each cell's momentum by its mass relative to total
                     pk = mom_per_cell*(d(ic+i,jc+j,kc+k)
     &                    + mass_per_cell)/m_nb/sqrt(real(dist, RKIND))
                     p_res(1) = p_res(1) + real(i, RKIND)*pk
                     p_res(2) = p_res(2) + real(j, RKIND)*pk
                     p_res(3) = p_res(3) + real(k, RKIND)*pk
                  endif
               endif
            enddo
         enddo
      enddo
c     Removing a bulk momentum drift would remove some kinetic energy
c     that should be deposited, so sum up how much and deposit it
c     as thermal energy later
      e_res = 0.5_RKIND*(p_res(1)**2 + p_res(2)**2 + p_res(3)**2)/m_res
c
c     Loop cells in particle-frame
c
      delta_mass =   mass_per_cell
      do k = -distrad, +distrad
         stepk = abs(k)
         do j = -distrad, +distrad
            stepj = stepk + abs(j)
            do i = -distrad, +distrad
               cellstep = stepj + abs(i)
               if (cellstep .le. diststep) then
               if (preSNmom .eq. 1_IKIND) then
c                 Calculate distance of this cell from host cell
                  dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
c                 Decompose the radial momentum onto Cartesian axes using
c                 the outward unit vector (i,j,k)/r where r = ||(i,j,k)||.
c                 At the center cell i=j=k=0 so delta_p components are zero;
c                 r=1 avoids division by zero in that case.
                  if (dist .eq. 0_IKIND) then
                     r = 1._RKIND
                  else
                     r = sqrt(real(dist, RKIND))
                  endif
c              Remove this cell's mass-weighted share of the bulk momentum
                  f_res = (d(ic+i,jc+j,kc+k) + delta_mass)/m_res
                  f_nb = (d(ic+i,jc+j,kc+k) + delta_mass)/m_nb
                  delta_pu = real(i, RKIND)/r * mom_per_cell * f_nb
     &                       - f_res*p_res(1)
                  delta_pv = real(j, RKIND)/r * mom_per_cell * f_nb
     &                       - f_res*p_res(2)
                  delta_pw = real(k, RKIND)/r * mom_per_cell * f_nb
     &                       - f_res*p_res(3)
                  delta_p = (delta_pu**2._RKIND + delta_pv**2._RKIND + 
     &               delta_pw**2._RKIND)**(0.5_RKIND)
c
c              delta_p's have units of mass*velocity/volume
c              delta_therm has units of mass*velocity^2/volume
#ifdef FORTRAN_DEBUG
               write(7,*) 'i,j,k, dist', i, j, k, dist
               write(7,*) 'p_before:',pu(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            pv(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            pw(ic+i,jc+j,kc+k)*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND
               write(7,*) 'delta_p:',delta_pu*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            delta_pv*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &            delta_pw*dx**3._RKIND*
     &            (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND
#endif
c
c     Add momentum
c           
                  pu_start = pu(ic+i,jc+j,kc+k)
                  pv_start = pv(ic+i,jc+j,kc+k)
                  pw_start = pw(ic+i,jc+j,kc+k)
                  p_start = (pu_start**2._RKIND + pv_start**2._RKIND +
     &               pw_start**2.)**(0.5_RKIND)
                  pu(ic+i,jc+j,kc+k) = 
     &                       pu_start + delta_pu
                  pv(ic+i,jc+j,kc+k) = 
     &                       pv_start + delta_pv
                  pw(ic+i,jc+j,kc+k) = 
     &                       pw_start + delta_pw
                  p_end = (pu(ic+i,jc+j,kc+k)**2._RKIND + 
     &                  pv(ic+i,jc+j,kc+k)**2._RKIND +
     &                  pw(ic+i,jc+j,kc+k)**2._RKIND)**(0.5_RKIND)
#ifdef FORTRAN_DEBUG
               write(7,*) 'v_start, v_final, mass:',
     &            p_start/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND,
     &            p_end/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND,
     &            (d(ic+i,jc+j,kc+k) + delta_mass)
     &              *dunits*(dx*x1)**3_IKIND/SolarMass
#endif
c              End momentum injection
               endif
c
c     Adjust thermal energy to account for new mass added
c
               dratio = d(i+ic,j+jc,k+kc)/
     &           (d(ic+i,jc+j,kc+k) + delta_mass)
               te(ic+i,jc+j,kc+k) = te(i+ic,j+jc,k+kc)*dratio
               if (idual .eq. 1_IKIND)
     &           ge(ic+i,jc+j,kc+k) = ge(ic+i,jc+j,kc+k)*dratio
c     
c     Metal feedback
c     metal field is metal fraction here, need to convert to density
c     before adding metal mass
c     
               if (imetal .eq. 1_IKIND) then
                  metal(ic+i,jc+j,kc+k) = 
     &               (metal(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &               metal_per_cell)/(d(ic+i,jc+j,kc+k) + delta_mass)
c                 Ensure the metallicity doesn't get too high
                  if (metal(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &              metal(ic+i,jc+j,kc+k) = 0.95_RKIND
c                 Dilute other metal fields by the mass added
                  if (itracksrc .eq. 1) then
                     metalSNII(ic+i,jc+j,kc+k) = 
     &                 (metalSNII(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k))
     &                 / (d(ic+i,jc+j,kc+k) + delta_mass)
                     metalSNIa(ic+i,jc+j,kc+k) = 
     &                 (metalSNIa(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k))
     &                 / (d(ic+i,jc+j,kc+k) + delta_mass)
                     if (useAGB .eq. 1_IKIND) metalAGB(ic+i,jc+j,kc+k) = 
     &                 (metalAGB(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k))
     &                 / (d(ic+i,jc+j,kc+k) + delta_mass)
                  endif
               endif
c
c     Add mass
c
               d(ic+i,jc+j,kc+k) = d(ic+i,jc+j,kc+k)
     &                       + delta_mass
c
c     End loop over cells in particle-frame
c
               endif
            enddo
         enddo
      enddo
      return
      end
c
c ==========================================================
c
c     Stochastically draw the number of SNe in this timestep
c
      subroutine draw_stochastic(nsn_ii, nsn_ia, dt)

      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      R_PREC    nsn_ii, nsn_ia, dt
c
c     Locals
c
      R_PREC    k, random
c
      k = aint(nsn_ii)
      call RANDOM_NUMBER(random)
      if (random .lt. nsn_ii - k) k = k + 1._RKIND
      nsn_ii = k
c
      k = aint(nsn_ia)
      call RANDOM_NUMBER(random)
      if (random .lt. nsn_ia - k) k = k + 1._RKIND
      nsn_ia = k
c
      return
      end