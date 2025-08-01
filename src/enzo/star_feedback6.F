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
     &               pSNtabMom)
c
c  RELEASES "STAR" PARTICLE MOMENTUM, MASS AND METALS
c
c  written by: Cassi Lochhaas, Claire Kopenhafer
c  date:      February 2025
c  	this is a copy of star_maker3mom.F that has been modifed to 
c       calculate the amount of momentum to inject from properties
c       of the host and surrounding cells, following Kimm & Cen (2014)
c       and Kimm et al. (2015).
c       The momentum is then injected into a cube of 27 cells centered
c       on the host cell of the particle.
c
c       Currently, the only momentum injection is due to Type II and
c       Type Ia SNe. Momentum injection from stellar winds is WIP.
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
c    use_tab_fdbk - flag to use tabular feedback
c    minit     - initial particle mass
c    ergSNII   - energy in ergs per SNII event
c    ergSNIa   - energy in ergs per SNIa event
c    itracksrc - track SNII & SNIa metals separately?
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
      INTG_PREC preSN, preSNmom, pSNntabZ, pSNntabAge
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    tabZ(ntabZ), tabAge(ntabAge)
      R_PREC    tabMass(4, ntabAge, ntabZ), tabMetal(4, ntabAge, ntabZ)
      R_PREC    tabEvents(2, ntabAge, ntabZ)
      R_PREC    pSNtabZ(pSNntabZ), pSNtabAge(pSNntabAge)
      R_PREC    pSNtabMass(pSNntabAge, pSNntabZ)
      R_PREC    pSNtabMetal(pSNntabAge, pSNntabZ)
      R_PREC    pSNtabMom(pSNntabAge, pSNntabZ)
c
c  Locals
c
      INTG_PREC n, ic, jc, kc, ip, jp, kp,
     &     i, j, k
      R_PREC mform, energy, energy_ii, energy_ia,
     &     m_eject, meject_ii, meject_ia, yield, minitial, xv1, xv2,
     &     distrad,distcells, mass_per_cell,
     &     mass_before, mass_dummy,
     &     kin_energy_before, mom_dummy, mom_inj, energy_inj,
     &     delta_ke, ke_injected, ke_after, kin_energy_inj,
     &     nsn, nsn_ii, nsn_ia, nsn_sto, nsn_ii_sto, nsn_ia_sto,
     &     mass_end, kin_energy_end, mom_start, mom_end, 
     &     kin_energy_dummy
      R_PREC xfcshift, yfcshift, zfcshift, xfc, yfc, zfc,
     &     xpos, ypos, zpos
      R_PREC ke_before(3,3,3), mom_per_cell(3,3,3)
      R_PREC energy_per_cell(3,3,3), energy_per_cell1(3,3,3)
      R_PREC Z_avg, n_avg, Z_floor,
     &	  realmass, mu_cell, fbuff, add_therm
      R_PREC n_cell, Z_cell, ke
      R_PREC te_sum, ke_sum, ge_sum, mom_sum, mass_sum, reduct_ratio,
     &     mass_term, vel_mag, delta_ge
      R_PREC mz_eject, mzeject_ii, mzeject_ia
      R_PREC metal_per_cell, metal_ii_per_cell, metal_ia_per_cell
      R_PREC u1(3,3,3), v1(3,3,3), w1(3,3,3), d1(3,3,3),
     &     te1(3,3,3), ge1(3,3,3), metal1(3,3,3)
      R_PREC nsn_grid(nx,ny,nz), meject_grid(nx,ny,nz),
     &     mzeject_grid(nx,ny,nz), mzeject_ii_grid(nx,ny,nz),
     &     mzeject_ia_grid(nx,ny,nz), up_grid(nx,ny,nz),
     &     vp_grid(nx,ny,nz), wp_grid(nx,ny,nz)
      R_PREC pSNmass_grid(nx,ny,nz), pSNmetal_grid(nx,ny,nz),
     &      pSNmom_grid(nx,ny,nz), pSNup_grid(nx,ny,nz),
     &      pSNvp_grid(nx,ny,nz), pSNwp_grid(nx,ny,nz)
      R_PREC up_avg, vp_avg, wp_avg, up_avg_w, vp_avg_w, wp_avg_w
      R_PREC wind_mass, wind_metal, wind_mom
      R_PREC ke_b_wind(3,3,3)
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
c     Ensure supernova energies are 1e51
c
      if (ergSNII .ne. 1e51_RKIND) then
         write(6,*) "Mechanical feedback (StarParticleFeedback = 64)" //
     &      " does not support supernova energies other than" //
     &      " 1e51 ergs."
         write(6,*) "If you want to use this feedback, set " //
     &      "StarFeedbackTabularSNIIEnergy = 1.e51 (its default value)."
         ERROR_MESSAGE
      endif
      if (ergSNIa .ne. 1e51_RKIND) then
         write(6,*) "Mechanical feedback (StarParticleFeedback = 64)" //
     &      " does not support supernova energies other than" //
     &      " 1e51 ergs."
         write(6,*) "If you want to use this feedback, set " //
     &      "StarFeedbackTabularSNIaEnergy = 1.e51 (its default value)."
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
            enddo
         enddo
      enddo
c
c     Assuming 3x3x3 cube, calculate cell distribution
c     NOTE: This feedback scheme is not written to be used with
c     any other cell distribution!!
c
      distrad = 3_IKIND
      distcells = distrad**3_IKIND
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
               write(7,*) 'warning: star particle out of grid',
     &              xp(n),yp(n),zp(n), xstart, ystart, zstart
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
               call pSN_mass(wind_mass, minit(n), metalf(n),
     &                          age, dt, t1, pSNtabMass, pSNtabZ,
     &                          pSNtabAge, pSNntabZ, pSNntabAge)
               call pSN_metal(wind_metal, minit(n), metalf(n),
     &                          age, dt, t1, pSNtabMetal, pSNtabZ,
     &                          pSNtabAge, pSNntabZ, pSNntabAge)
               if (preSNmom .eq. 1_IKIND) then
                  call pSN_mom(wind_mom, minit(n), metalf(n), age,
     &                          dt, t1, vunits, pSNtabMom, pSNtabZ,
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
c     If doing stochastic feedback, use the expected number of SNe in
c     this timestep to draw from a Poisson distribution to get an integer
c     number
c
               if (stochastic .eq. 1_IKIND) then
                  nsn_timestep = 1._RKIND
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
c     with 0.25*current mass to ensure no more than 25% of the particle's mass
c     can be ejected at once.
c
                  if (nsn_ii_sto .gt. 0._RKIND) then
                     meject_ii = min((meject_ii / nsn_ii * nsn_ii_sto)* 
     &                           mp(n)/minit(n), 0.25_RKIND*mp(n))
                  else
                     meject_ii = 0._RKIND
                  endif
                  if (nsn_ia_sto .gt. 0._RKIND) then
                     meject_ia = min((meject_ia / nsn_ia * nsn_ia_sto)*
     &                           mp(n)/minit(n), 0.25_RKIND*mp(n))
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
c     Metals done similarly to mass, but ensure no more than 2% of current
c     mass ejected as metals
c
                  if (imetal .eq. 1_IKIND) then
                     if (nsn_ii_sto .gt. 0._RKIND) then
                        mzeject_ii = min((mzeject_ii / nsn_ii * 
     &                    nsn_ii_sto)*mp(n)/minit(n), 0.02_RKIND*mp(n))
                     else
                        mzeject_ii = 0._RKIND
                     endif
                     if (nsn_ia_sto .gt. 0._RKIND) then
                        mzeject_ia = min((mzeject_ia / nsn_ia * 
     &                    nsn_ia_sto)*mp(n)/minit(n), 0.02_RKIND*mp(n))
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
               endif ! end stochastic SNe
c
               energy = nsn * 10._RKIND**51._RKIND / (dunits*
     &                  x1**3._RKIND*vunits**2._RKIND*dx**3._RKIND)
#ifdef FORTRAN_DEBUG
               write(7,*) 'nsn,m_eject,energy:',nsn,m_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              energy*(dunits*x1**3._RKIND*vunits**2._RKIND*
     &               dx**3._RKIND)
#endif
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
               if (preSN .eq. 1_IKIND) mp(n) = mp(n) - wind_mass
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
               energy = nsn * 10._RKIND**51._RKIND / (dunits*
     &                  x1**3._RKIND*vunits**2._RKIND*dx**3._RKIND)
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
               write(7,*) 'nsn,m_eject,energy:',nsn,m_eject*
     &              dunits*x1**3._RKIND/SolarMass*dx**3_IKIND,
     &              energy*(dunits*x1**3._RKIND*vunits**2._RKIND*
     &               dx**3._RKIND)
#endif
            endif ! end Cen & Ostriker method
c
c	   Set center of feedback zone
c
            xfc = xp(n)
            yfc = yp(n)
            zfc = zp(n)
            fbuff = ibuff + 2._RKIND
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
c     Put nsn, meject, mzeject, and weighted particle velocity values on grid at this location
c     Weight average velocity for SNe by number of SNe going off in each particle
c
            nsn_grid(ic,jc,kc) = nsn_grid(ic,jc,kc) + nsn
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
c     Loop over the feedback grids, avoiding the edges
c
      do kc=2_IKIND, nz-1_IKIND
         do jc=2_IKIND, ny-1_IKIND
            do ic=2_IKIND, nx-1_IKIND
c
c ====================================================================
c     BEGIN PRE-SUPERNOVA FEEDBACK
c
c     If injecting momentum from pre-SN feedback:
c
               if (preSN .eq. 1_IKIND .and. 
     &               pSNmass_grid(ic,jc,kc) .gt. 1e-15_RKIND) then
                  if (preSNmom .eq. 1_IKIND) then
c
c     Compute the kinetic energy in the affected region
c     before momentum is added.  This
c     is needed at the end of the calculation to update the
c     total energy (te) field.
c
               do k = -1_IKIND, +1_IKIND
                  do j = -1_IKIND, +1_IKIND
                     do i = -1_IKIND, +1_IKIND
                        ke_b_wind(i+2, j+2, k+2) = 
     &                      0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                              (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                               v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                               w(ic+i ,jc+j ,kc+k)**2_IKIND)
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
     &                    +1_IKIND)
c     Sum mass, energy, and momentum before
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc,
     &            mass_b_wind, kin_energy_b_wind, mom_b_wind)
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
               do k = -1_IKIND, 1_IKIND
                  do j = -1_IKIND, 1_IKIND
                     do i = -1_IKIND, 1_IKIND
                        u1(i+2,j+2,k+2) = 0._RKIND
                        v1(i+2,j+2,k+2) = 0._RKIND
                        w1(i+2,j+2,k+2) = 0._RKIND
                        d1(i+2,j+2,k+2) = d(i+ic,j+jc,k+kc)
                        ge1(i+2,j+2,k+2) = 0._RKIND
                        te1(i+2,j+2,k+2) = 0._RKIND
                        metal1(i+2,j+2,k+2) = 0._RKIND
                     enddo
                  enddo
               enddo
c
c     Now add mass and momentum terms to local dummy fields
c              
               mass_per_cell_w = pSNmass_grid(ic,jc,kc)/distcells
               mom_per_cell_w = pSNmom_grid(ic,jc,kc)
     &                              /(distcells-1_IKIND)
               metal_per_cell_w = pSNmetal_grid(ic,jc,kc)/distcells
               call add_preSN_feedback(u1, v1, w1, d1, ge1, te1, metal1, 
     &                        3_IKIND, 3_IKIND, 3_IKIND,
     &                        2_IKIND, 2_IKIND, 2_IKIND,
     &                        imetal, imulti_metals, idual,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell_w, mom_per_cell_w,
     &                        metal_per_cell_w, dx, x1, dunits, vunits)
c
c     Calculate mass and energy in dummy fields
c
               call sum_energy_mom(u1, v1, w1, d1,
     &                        3_IKIND, 3_IKIND, 3_IKIND,
     &                        2_IKIND, 2_IKIND, 2_IKIND,
     &                        mass_dummy_w, kin_energy_dummy_w, 
     &                        mom_dummy_w)
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
     &                        imetal, imulti_metals, idual,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell_w, mom_per_cell_w,
     &                        metal_per_cell_w,
     &                        dx, x1, dunits, vunits)
c
c     Calculate the mass and energy in the affected region again.
c
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &                         ic, jc, kc,
     &                         mass_end_w, kin_energy_end_w, mom_end_w)
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
c     the remaining energy as thermal into the host cell.
c 
               add_therm_w = 0._RKIND
               kin_energy_inj_w = kin_energy_end_w - kin_energy_b_wind
               if (kin_energy_inj_w .lt. kin_energy_dummy_w) then
                  add_therm_w = (kin_energy_dummy_w-kin_energy_inj_w)
                  te(ic,jc,kc) = te(ic,jc,kc) + 
     &                           add_therm_w/d(ic,jc,kc)
                  if (idual .eq. 1_IKIND) ge(ic,jc,kc) = ge(ic,jc,kc) + 
     &                           add_therm_w/d(ic,jc,kc)
#ifdef FORTRAN_DEBUG
                  write(7,*) "canceled wind mom injected as thermal:",
     &                  add_therm_w*dunits*
     &                  vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
#endif
               endif

c
c           Convert momenta back to velocities and transform back to lab frame
c
               call momentum_convert(u, v, w, d, up_avg_w, vp_avg_w, 
     &               wp_avg_w, nx, ny, nz, ic, jc, kc, -1_IKIND)
c
c           Add the change in the kinetic energy to the total energy field
c     
               do k = -1_IKIND, +1_IKIND
                  do j = -1_IKIND, +1_IKIND
                     do i = -1_IKIND, +1_IKIND
                        vel_mag = (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                    v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                    w(ic+i ,jc+j ,kc+k)**2_IKIND)**0.5_RKIND
                        ke_after_w = 0.5_RKIND * d(ic+i, jc+j, kc+k) * 
     &                    vel_mag**2_IKIND
                        delta_ke_w = ke_after_w - 
     &                        ke_b_wind(i+2,j+2,k+2)
                        te(ic+i,jc+j,kc+k) = te(ic+i,jc+j,kc+k) +
     &                     delta_ke_w/d(ic+i,jc+j,kc+k)
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
     &            mass_b_wind, kin_energy_b_wind, mom_b_wind)
               write(7,*) "m before wind:", mass_b_wind*
     &               (dunits*(x1*dx)**3._RKIND)/SolarMass
c
#endif
               mass_per_cell_w = pSNmass_grid(ic,jc,kc)/distcells
               metal_per_cell_w = pSNmetal_grid(ic,jc,kc)/distcells
               call add_preSN_feedback(u, v, w, d, ge, te, metal,
     &                         nx, ny, nz, ic, jc, kc,
     &                        imetal, imulti_metals, idual,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell_w, 0._RKIND,
     &                        metal_per_cell_w,
     &                        dx, x1, dunits, vunits)
#ifdef FORTRAN_DEBUG
c     Sum mass, energy, and momentum after. Note mom and KE calculations
c     will be wrong here because velocities haven't been converted to
c     momentum, but we only care about mass for this debugging check
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc,
     &            mass_end_w, kin_energy_end_w, mom_end_w)
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
c     Compute average density and metallicity around this cell
c     
               Z_avg = 0._RKIND
               n_avg = 0._RKIND
               do k = -1_IKIND,1_IKIND
                  do j = -1_IKIND,1_IKIND
                     do i = -1_IKIND,1_IKIND
                        Z_cell = metal(ic+i,jc+j,kc+k)
                        mu_cell = mu(ic+i,jc+j,kc+k)
c
c     The following adjustments to Z_cell and mu_cell are band-aids to
c     handle situations when the chemistry gets weird and the metallicity
c     is unphysically high or mu as calculated from the baryon fields in
c     Grid_StarParticleHandler.C is somehow negative, NaN, or outside of the
c     expected range from ~0.5 to ~16. The below calculation
c     of mu assumes primordial relative abundance between H and He (He is
c     1/3 of H), that the gas is completely neutral, and that all metals
c     are oxygen atoms, for simplicity. Note that what is called "mu" here
c     and in the mu_field passed to this function are actually the inverse
c     mean molecular weight!!!
c
                        if (Z_cell .gt. 0.95_RKIND) 
     &                    Z_cell = 0.95_RKIND
                        if ((mu_cell .lt. 0.0625_RKIND) .or. 
     &                    (mu_cell .gt. 2._RKIND) .or.
     &                    (mu_cell .ne. mu_cell))
     &                    mu_cell = 0.8125_RKIND*(1._RKIND-Z_cell) + 
     &                     Z_cell/16._RKIND
                        Z_avg = Z_avg + Z_cell/0.02_RKIND
                        n_avg = n_avg + 
     &                     d(ic+i,jc+j,kc+k)*dunits*
     &                     mu_cell/mass_h
#ifdef FORTRAN_DEBUG
                        write(7,*) 'mu, Z, d, n',
     &                       mu_cell, Z_cell, d(ic+i,jc+j,kc+k)*dunits,
     &                       d(ic+i,jc+j,kc+k)*dunits*mu_cell/mass_h
#endif
                     enddo
                  enddo
               enddo
c               
c     The Kimm & Cen (2014) scheme (eq A5) requires a floor on Z_avg
c
               Z_floor = max(Z_avg/distcells, 0.01_RKIND)
               n_avg = n_avg/distcells
#ifdef FORTRAN_DEBUG
               write(7,*) 'Z_avg, n_avg', Z_floor, n_avg
#endif
c
c     Compute the mass and momentum that this cell is responsible for
c     injecting
c
               mass_per_cell = meject_grid(ic,jc,kc)/distcells
               call get_mom_per_cell(d, mass_per_cell,
     &                  distcells, mom_per_cell, energy_per_cell,
     &                  nsn_grid(ic,jc,kc), d(ic,jc,kc), dx, x1,
     &                  n_avg, Z_floor, dunits, vunits, t1,
     &                  nx, ny, nz, ic, jc, kc,
     &                  mom_mult)
c           mom_per_cell has units of mass*velocity/volume
c
c     If not injecting any momentum or energy, skip
c
               mom_inj = 0._RKIND
               energy_inj = 0._RKIND
               do k = 1_IKIND,3_IKIND
                  do j = 1_IKIND,3_IKIND
                     do i = 1_IKIND,3_IKIND
                        mom_inj = mom_inj + mom_per_cell(i, j, k)
                        energy_inj = energy_inj + 
     &                      energy_per_cell(i, j, k)
                     enddo
                  enddo
               enddo
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
               do k = -1_IKIND, +1_IKIND
                  do j = -1_IKIND, +1_IKIND
                     do i = -1_IKIND, +1_IKIND
                        ke_before(i+2, j+2, k+2) = 
     &                      0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                              (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                               v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                               w(ic+i ,jc+j ,kc+k)**2_IKIND)
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
     &                    +1_IKIND)
c
c           u,v,w have units of mass*velocity/volume
c           Now u,v,w are the explosion-relative *momentum* directions
c
c     Sum mass, energy, and momentum before
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &            ic, jc, kc,
     &            mass_before, kin_energy_before, mom_start)
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
               do k = -1_IKIND, 1_IKIND
                  do j = -1_IKIND, 1_IKIND
                     do i = -1_IKIND, 1_IKIND
                        u1(i+2,j+2,k+2) = 0._RKIND
                        v1(i+2,j+2,k+2) = 0._RKIND
                        w1(i+2,j+2,k+2) = 0._RKIND
                        d1(i+2,j+2,k+2) = d(i+ic,j+jc,k+kc)
                        ge1(i+2,j+2,k+2) = 0._RKIND
                        te1(i+2,j+2,k+2) = 0._RKIND
                        metal1(i+2,j+2,k+2) = 0._RKIND
                        energy_per_cell1(i+2,j+2,k+2) = 0._RKIND
                     enddo
                  enddo
               enddo
c
c     Now add mass and momentum terms to local dummy fields
c     Metals and thermal energy zeroed out because we don't care about them in the dummy
c
               call add_feedback_mech(u1, v1, w1, d1, ge1, te1, metal1, 
     &                        3_IKIND, 3_IKIND, 3_IKIND,
     &                        2_IKIND, 2_IKIND, 2_IKIND,
     &                        0_IKIND, imulti_metals, idual,
     &                        itracksrc, kick_cap,
     &                        0._RKIND, 0._RKIND, 0._RKIND,
     &                        metalSNII, metalSNIa,
     &                        mass_per_cell, mom_per_cell,
     &                        energy_per_cell1,
     &                        dx, x1, dunits, vunits)
c
c     Calculate mass and energy in dummy fields
c
               call sum_energy_mom(u1, v1, w1, d1,
     &                        3_IKIND, 3_IKIND, 3_IKIND,
     &                        2_IKIND, 2_IKIND, 2_IKIND,
     &                        mass_dummy, kin_energy_dummy, mom_dummy)
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
     &                        itracksrc, kick_cap,
     &                        metal_per_cell, 
     &                        metal_ia_per_cell, 
     &                        metal_ii_per_cell,
     &                        metalSNII, metalSNIa,
     &                        mass_per_cell, mom_per_cell,
     &                        energy_per_cell,
     &                        dx, x1, dunits, vunits)
c
c     Calculate the mass and energy in the affected region again.
c
               call sum_energy_mom(u, v, w, d, nx, ny, nz,
     &                          ic, jc, kc,
     &                          mass_end, kin_energy_end, mom_end)
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
c     (because of momentum cancellation due to non-static gas), add
c     the remaining energy as thermal into the host cell.
c
               kin_energy_inj = kin_energy_end - kin_energy_before
               if (kin_energy_inj .lt. kin_energy_dummy) then
                  add_therm = (kin_energy_dummy-kin_energy_inj)
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
c           Convert momenta back to velocities and transform back to lab frame
c
               call momentum_convert(u, v, w, d, up_avg, vp_avg, 
     &               wp_avg, nx, ny, nz, ic, jc, kc, -1_IKIND)
c
c           Add the change in the kinetic energy to the total energy field
c     
               ke_injected = 0._RKIND
               te_sum = 0._RKIND
               ge_sum = 0._RKIND
               ke_sum = 0._RKIND
               mom_sum = 0._RKIND
               mass_sum = 0._RKIND
               do k = -1_IKIND, +1_IKIND
                  do j = -1_IKIND, +1_IKIND
                     do i = -1_IKIND, +1_IKIND
                        vel_mag = (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                    v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                    w(ic+i ,jc+j ,kc+k)**2_IKIND)**0.5_RKIND
                        ke_after = 0.5_RKIND * d(ic+i, jc+j, kc+k) * 
     &                    vel_mag**2_IKIND
                        delta_ke = ke_after - ke_before(i+2,j+2,k+2)
                        te(ic+i,jc+j,kc+k) = te(ic+i,jc+j,kc+k) +
     &                     delta_ke/d(ic+i,jc+j,kc+k)
                        ke_injected = ke_injected + delta_ke
                        mass_term = d(ic+i,jc+j,kc+k)
                        ke_sum = ke_sum + ke_after
                        mom_sum = mom_sum + mass_term*vel_mag
                        mass_sum = mass_sum + mass_term
c                       ke_after has units of mass*velocity^2/volume
c                       te has units of mass*velocity^2/volume/mass*volume = velocity^2
c                       te is *specific* energy (energy/mass)
                     enddo
                  enddo
               enddo
c
c          Inject the energy lost due to colliding flows as thermal energy
c          spread through the injection region.
c
               delta_ge = add_therm/distcells
               do k = -1_IKIND, +1_IKIND
                  do j = -1_IKIND, +1_IKIND
                     do i = -1_IKIND, +1_IKIND
                        te(ic+i,jc+j,kc+k) = te(ic+i,jc+j,kc+k) + 
     &                    delta_ge/d(ic+i,jc+j,kc+k)
                        ke = 0.5_RKIND*
     &                    (u(ic+i,jc+j,kc+k)**2_IKIND +
     &                    v(ic+i,jc+j,kc+k)**2_IKIND +
     &                    w(ic+i,jc+j,kc+k)**2_IKIND)
                        if (idual .eq. 1_IKIND) then
                           ge(ic+i,jc+j,kc+k) = ge(ic+i,jc+j,kc+k) + 
     &                       delta_ge/d(ic+i,jc+j,kc+k)
#ifdef FORTRAN_DEBUG
                           write(7,*) 'ge, te, ke after injection:',
     &                       ge(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k)*
     &                       dunits*vunits**2._RKIND*(x1*dx)**3._RKIND,
     &                       te(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k)*
     &                       dunits*vunits**2._RKIND*(x1*dx)**3._RKIND,
     &                       ke*d(ic+i,jc+j,kc+k)*
     &                       dunits*vunits**2._RKIND*(x1*dx)**3._RKIND
#endif
                           ge_sum = ge_sum +
     &                       ge(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k)
                        endif
                        te_sum = te_sum + 
     &                    te(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k)
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
     &             add_therm*dunits*
     &             vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND,
     &             pSNmass_grid(ic,jc,kc)*
     &             (dunits*(x1*dx)**3._RKIND)/SolarMass,
     &             pSNmetal_grid(ic,jc,kc)*(dunits*(x1*dx)**3._RKIND)/
     &             SolarMass,
     &             pSNmom_grid*(dunits*(dx*x1)**3._RKIND)*
     &             vunits/SolarMass/1e5_RKIND,
     &             add_therm_w*dunits*
     &             vunits**2._RKIND*dx**3._RKIND*x1**3._RKIND
               endif
c
 100           continue
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
c ==========================================================
c
c     Calculate the amount of momentum to add to each neighboring cell
c     Formulae taken from Kimm & Cen (2014), Appendix A
c
      subroutine get_mom_per_cell(d, mass_per_cell, 
     &               distcells, mom_per_cell, energy_per_cell,
     &               nsn, d_fbck_cell, dx, x1,
     &               n_avg, Z_floor, dunits, vunits, t1,
     &               nx, ny, nz, ic, jc, kc,
     &               mom_mult)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      R_PREC    d(nx,ny,nz), dx
      R_PREC    t1, Z_floor, n_avg
      R_PREC    mom_per_cell(3_IKIND, 3_IKIND, 3_IKIND)
      R_PREC    energy_per_cell(3_IKIND, 3_IKIND, 3_IKIND)
      R_PREC    mass_per_cell, distcells
      R_PREC    nsn, x1, dunits, d_fbck_cell, vunits
      R_PREC    mom_mult
c
c     Locals
c
      INTG_PREC i, j, k, dist
      R_PREC    dMswept, chi, chi_th, fe, pSN, ergs_51_sqr, ke_from_p
      R_PREC    den_cell, Mej, dMej, Nnbors, betaSN, SN_energy
c
        ergs_51_sqr = 10._RKIND**(51._RKIND/2._RKIND)
        Nnbors = distcells - 1._RKIND ! 26
        betaSN = 1._RKIND/distcells   ! 1/27
        Mej = mass_per_cell*distcells*dunits*(x1*dx)**3._RKIND
        SN_energy = nsn*10**51._RKIND/Nnbors
        chi_th = 69.58_RKIND*(nsn)**(-2._RKIND/17._RKIND) * 
     &               (n_avg)**(-4._RKIND/17._RKIND) * 
     &               (Z_floor)**(-0.28_RKIND)
#ifdef FORTRAN_DEBUG
        write(7,*) 'chi_th:',chi_th
#endif
        do k = -1_IKIND, +1_IKIND
            do j = -1_IKIND, +1_IKIND
               do i = -1_IKIND, +1_IKIND
c
c       Compute the mass fraction of swept-up to ejected mass.
c
                  den_cell = d(ic+i,jc+j,kc+k)*dunits
                  dMej = (1._RKIND-betaSN)*Mej/Nnbors
                  dMswept = den_cell * 
     &                     (dx*x1/2._RKIND)**(3._RKIND) + 
     &                     ((1._RKIND-betaSN)*
     &                     (d_fbck_cell*dunits)*(dx*x1)**3._RKIND) /
     &                     Nnbors + dMej
c
                  chi = dMswept/dMej
#ifdef FORTRAN_DEBUG
                  write(7,*) 'dMej,dMswept,chi:',dMej,dMswept,chi
#endif
c
c       Compare chi to threshold and set momentum injection.
c
                  if (chi .lt. chi_th) then
                     fe = 1._RKIND - ((chi - 1._RKIND)/
     &                        (3._RKIND*(chi_th - 1._RKIND)))
                     pSN = mom_mult * (2._RKIND*chi * Mej * fe *
     &                        nsn)**(0.5_RKIND) * ergs_51_sqr
c
c        This is the resolved case, so need to inject thermal energy
c        to get up to 10^51 per SNe. Use pSN to calculate kinetic energy
c        and inject thermal into this cell equal to total - kinetic
                     ke_from_p = 0.5_RKIND*(pSN/Nnbors)**2. / 
     &                      (den_cell*(dx*x1)**3._RKIND + dMej)
                     if (ke_from_p .lt. SN_energy) then
                           energy_per_cell(i+2,j+2,k+2) = 
     &                     (SN_energy - ke_from_p)
     &                     /(dunits*(dx*x1)**3._RKIND*vunits**2._RKIND)
                     else
                        energy_per_cell(i+2,j+2,k+2) = 0._RKIND
                     endif
#ifdef FORTRAN_DEBUG
                     write(7,*) 'lt chi_th'
                     write(7,*) 'pSN,ke_from_p,SN_energy,therm:',
     &                    pSN/SolarMass/1e5_RKIND, 
     &                    ke_from_p,
     &                    SN_energy,
     &                    energy_per_cell(i+2,j+2,k+2)*
     &                    (dunits*(dx*x1)**3._RKIND*vunits**2._RKIND)
#endif
                  else
                     pSN = mom_mult * 3.e10_RKIND*SolarMass * 
     &                        (nsn)**(16._RKIND/17._RKIND) * 
     &                        (n_avg)**(-2._RKIND/17._RKIND) * 
     &                        (Z_floor)**(-0.14_RKIND)
                     energy_per_cell(i+2,j+2,k+2) = 0._RKIND
#ifdef FORTRAN_DEBUG
                     write(7,*) 'gt chi_th, pSN:',
     &                        pSN/SolarMass/1e5_RKIND
#endif
                  endif
c              pSN has units of cgs-mass*cgs-velocity, want it to have units of code-mass*code-velocity/code-volume
c              conversion of mass from code units is (dunits*(dx*x1)**3._RKIND) (density times volume)
                     pSN = pSN/(dunits*(dx*x1)**3._RKIND)/vunits/Nnbors
c              dunits = cgs-density/code-density = cgs-mass*code-volume/(cgs-volume*code-mass)
c              x1^3 = cgs-volume/code-volume
c              vunits = cgs-velocity/code-velocity
c              dx^3 = code-volume
c              dunits*(dx*x1)**3._RKIND = cgs-mass*code-volume/(cgs-volume*code-mass)*code-volume*cgs-volume/code-volume
c              dunits*(dx*x1)**3._RKIND = cgs-mass*code-volume/code-mass
c              pSN = cgs-mass*cgs-velocity/(cgs-mass*code-volume/code-mass)/(cgs-velocity/code-velocity)
c              pSN = code-mass*code-velocity/code-volume
c              pSN now has units code-mass*code-velocity/code-volume
c              So pSN is actually in the right units after all!
c
c        If cell is in center of cube, no momentum or energy.
c
                     dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
                     if (dist .eq. 0_IKIND) then
                        pSN = 0._RKIND
                        energy_per_cell(i+2,j+2,k+2) = 0._RKIND
                     endif
                  mom_per_cell(i+2,j+2,k+2) = pSN
c
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
     &                    nx, ny, nz, ic, jc, kc, idir)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      INTG_PREC idir
      R_PREC    d(nx, ny, nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    up, vp, wp
c
c     Locals
c
      INTG_PREC i, j, k
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
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
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
     &                           mass_sum, kin_energy_sum, mom_sum)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    d(nx,ny,nz)
      R_PREC    mass_sum, kin_energy_sum
      R_PREC    mom_sum
c
c     Locals
c
      INTG_PREC i, j, k
      R_PREC    mass_term, mom_term, kin_energy
c
      mass_sum = 0._RKIND
      kin_energy_sum = 0._RKIND
      mom_sum = 0._RKIND
c
c     Sum mass and energy      
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
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
     &                        itracksrc, kick_cap,
     &                        mzeject, mzeject_ia, mzeject_ii,
     &                        metalSNII, metalSNIa,
     &                        mass_per_cell, mom_per_cell,
     &                        energy_per_cell,
     &                        dx, x1, dunits, vunits)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual, dist
      INTG_PREC imetal, imulti_metals, kick_cap
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_per_cell, mom_per_cell(3_IKIND,3_IKIND,3_IKIND)
      R_PREC    energy_per_cell(3_IKIND,3_IKIND,3_IKIND)
      R_PREC    mzeject, mzeject_ia, mzeject_ii
      R_PREC    metalSNII(nx,ny,nz), metalSNIa(nx,ny,nz)
      R_PREC    dunits, x1, vunits, dx, mult
      R_PREC    p_start, p_end
      INTG_PREC    itracksrc
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw,
     &          delta_p, vel_kick, dratio
      R_PREC    pu_start, pu_deposit, pv_start, pv_deposit, pw_start,
     &          pw_deposit
c
c     Loop cells in particle-frame
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c              Calculate distance of this cell from host cell
               dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
c              Multiply the momentum deposited in corners or edges of
c              cubes by the correct geometric factor for each x,y,z component
               if (dist .eq. 0_IKIND) then
                  mult = 0._RKIND
               else if (dist .eq. 1_IKIND) then
                  mult = 1._RKIND
               else if (dist .eq. 2_IKIND) then
                  mult = 1._RKIND/(2._RKIND**0.5_RKIND)
               else if (dist .eq. 3_IKIND) then
                  mult = 1._RKIND/(3._RKIND**0.5_RKIND)
               endif
               delta_mass =   mass_per_cell
c              delta_mass has units of mass/volume
               delta_pu   =  i*mult*mom_per_cell(i+2,j+2,k+2)
               delta_pv   =  j*mult*mom_per_cell(i+2,j+2,k+2)
               delta_pw   =  k*mult*mom_per_cell(i+2,j+2,k+2)
               delta_p = (delta_pu**2._RKIND + delta_pv**2._RKIND + 
     &            delta_pw**2._RKIND)**(0.5_RKIND)
c              If delta_p leads to a velocity kick of more than
c              1000 km/s, limit to 1000
               vel_kick = delta_p/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND
               if (kick_cap .eq. 1_IKIND) then
                  if (vel_kick .gt. 1000._RKIND) then
                     delta_p = 1000._RKIND/(vunits/1e5_RKIND)*
     &                   (d(ic+i,jc+j,kc+k) + delta_mass)
                     delta_pu = i*mult*delta_p
                     delta_pv = j*mult*delta_p
                     delta_pw = k*mult*delta_p
                  endif
               endif
#ifdef FORTRAN_DEBUG
                  write(7,*) 'vel_kick, delta_p, delta_v:',
     &              vel_kick,
     &              delta_p*dx**3._RKIND*
     &              (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &              delta_p/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND
#endif
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
c     Metal feedback (note that in this function gas metal is
c     a fraction (rho_metal/rho_gas)
c     
               if (imetal .eq. 1_IKIND) then
                  metal(ic+i,jc+j,kc+k) = 
     &               (metal(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &               mzeject) /(d(ic+i, jc+j, kc+k) + delta_mass)
c                 Ensure the metallicity doesn't get too high
                  if (metal(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &              metal(ic+i,jc+j,kc+k) = 0.95_RKIND
                  if (itracksrc .eq. 1) then
                     metalSNII(ic+i,jc+j,kc+k) =
     &                 (metalSNII(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &                 mzeject_ii) /(d(ic+i, jc+j, kc+k) + delta_mass)
c                 Put same cap on other metal fields
                     if (metalSNII(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &                 metalSNII(ic+i,jc+j,kc+k) = 0.95_RKIND
                     metalSNIa(ic+i,jc+j,kc+k) =
     &                 (metalSNIa(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &                 mzeject_ia) /(d(ic+i, jc+j, kc+k) + delta_mass)
                     if (metalSNIa(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &                 metalSNIa(ic+i,jc+j,kc+k) = 0.95_RKIND
                  endif
               endif
c
c     Add mass
c
               d(ic+i,jc+j,kc+k) = d(ic+i,jc+j,kc+k)
     &                       + delta_mass
c
c     Add thermal energy for any cells where Sedov-Taylor phase is resolved
c     (calculated in get_mom_per_cell)
c
               te(ic+i,jc+j,kc+k) = te(i+ic,j+jc,k+kc) + 
     &               energy_per_cell(i+2,j+2,k+2)/d(ic+i,jc+j,kc+k)
               if (idual .eq. 1_IKIND)
     &           ge(ic+i,jc+j,kc+k) = ge(ic+i,jc+j,kc+k)+ 
     &               energy_per_cell(i+2,j+2,k+2)/d(ic+i,jc+j,kc+k)
c     
c     End loop over cells in particle-frame
c     
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
     &                        imetal, imulti_metals, idual,
     &                        kick_cap, preSNmom,
     &                        mass_per_cell, mom_per_cell,
     &                        metal_per_cell,
     &                        dx, x1, dunits, vunits)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc, idual, dist
      INTG_PREC imetal, imulti_metals, kick_cap, preSNmom
      R_PREC    d(nx, ny, nz), metal(nx,ny,nz), ge(nx,ny,nz)
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz), te(nx,ny,nz)
      R_PREC    mass_per_cell, mom_per_cell
      R_PREC    metal_per_cell
      R_PREC    dunits, x1, vunits, dx, mult
      R_PREC    p_start, p_end
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    delta_mass, delta_pu, delta_pv, delta_pw,
     &          delta_p, vel_kick, dratio
      R_PREC    pu_start, pu_deposit, pv_start, pv_deposit, pw_start,
     &          pw_deposit
c
c     Loop cells in particle-frame
c
      delta_mass =   mass_per_cell
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
               if (preSNmom .eq. 1_IKIND) then
c                 Calculate distance of this cell from host cell
                  dist = i**2_IKIND + j**2_IKIND + k**2_IKIND
c                 Multiply the momentum deposited in corners or edges of
c                 cubes by the correct geometric factor for each x,y,z component
                  if (dist .eq. 0_IKIND) then
                     mult = 0._RKIND
                  else if (dist .eq. 1_IKIND) then
                     mult = 1._RKIND
                  else if (dist .eq. 2_IKIND) then
                     mult = 1._RKIND/(2._RKIND**0.5_RKIND)
                  else if (dist .eq. 3_IKIND) then
                     mult = 1._RKIND/(3._RKIND**0.5_RKIND)
                  endif
c                 delta_mass has units of mass/volume
                  delta_pu   =  i*mult*mom_per_cell
                  delta_pv   =  j*mult*mom_per_cell
                  delta_pw   =  k*mult*mom_per_cell
                  delta_p = (delta_pu**2._RKIND + delta_pv**2._RKIND + 
     &               delta_pw**2._RKIND)**(0.5_RKIND)
c                 If delta_p leads to a velocity kick of more than
c                 1000 km/s, limit to 1000
                  vel_kick = delta_p/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND
                  if (kick_cap .eq. 1_IKIND) then
                     if (vel_kick .gt. 1000._RKIND) then
                        delta_p = 1000._RKIND/(vunits/1e5_RKIND)*
     &                   (d(ic+i,jc+j,kc+k) + delta_mass)
                        delta_pu = i*mult*delta_p
                        delta_pv = j*mult*delta_p
                        delta_pw = k*mult*delta_p
                     endif
                  endif
#ifdef FORTRAN_DEBUG
                  write(7,*) 'vel_kick, delta_p, delta_v:',
     &              vel_kick,
     &              delta_p*dx**3._RKIND*
     &              (dunits*x1**3._RKIND)*vunits/SolarMass/1e5_RKIND,
     &              delta_p/(d(ic+i,jc+j,kc+k) + delta_mass)*
     &                    vunits/1e5_RKIND
#endif
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
c     Metal feedback (note that in this function gas metal is
c     a fraction (rho_metal/rho_gas)
c     
               if (imetal .eq. 1_IKIND) then
                  metal(ic+i,jc+j,kc+k) = 
     &               (metal(ic+i,jc+j,kc+k)*d(ic+i,jc+j,kc+k) + 
     &               metal_per_cell) /(d(ic+i, jc+j, kc+k) + delta_mass)
c                 Ensure the metallicity doesn't get too high
                  if (metal(ic+i,jc+j,kc+k) .gt. 0.95_RKIND)
     &              metal(ic+i,jc+j,kc+k) = 0.95_RKIND
               endif
c
c     Add mass
c
               d(ic+i,jc+j,kc+k) = d(ic+i,jc+j,kc+k)
     &                       + delta_mass
c     
c     End loop over cells in particle-frame
c     
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
      R_PREC    k, P, L_ii, L_ia, random
c
      L_ii = EXP(-nsn_ii)
      L_ia = EXP(-nsn_ia)
      k = 0._RKIND
      P = 1._RKIND
c
      do
         call RANDOM_NUMBER(random)
         P = P * random
         if (P .lt. L_ii) goto 20
         k = k + 1._RKIND
      enddo
 20   continue
      nsn_ii = k
c
      P = 1._RKIND
      k = 0._RKIND
c
      do
         call RANDOM_NUMBER(random)
         P = P * random
         if (P .lt. L_ia) goto 30
         k = k + 1._RKIND
      enddo
 30   continue
      nsn_ia = k
c
      return
      end