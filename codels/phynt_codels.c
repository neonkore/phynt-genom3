/*
 * Copyright (c) 2018 LAAS/CNRS
 * All rights reserved.
 *
 * Redistribution and use  in source  and binary  forms,  with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *   1. Redistributions of  source  code must retain the  above copyright
 *      notice and this list of conditions.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice and  this list of  conditions in the  documentation and/or
 *      other materials provided with the distribution.
 *
 *                                      Anthony Mallet on Mon Jun 11 2018
 */
#include "acphynt.h"

#include <sys/time.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include "phynt_c_types.h"
#include "codels.h"


/* --- Attribute set_af_parameters -------------------------------------- */

/** Validation codel phynt_set_af_parameters of attribute set_af_parameters.
 *
 * Returns genom_ok.
 * Throws .
 */
genom_event
phynt_set_af_parameters(const double J[9], const genom_context self)
{
  phynt_admittance_J(J);
  return genom_ok;
}


/* --- Function set_state ----------------------------------------------- */

/** Codel phynt_set_state of function set_state.
 *
 * Returns genom_ok.
 */
genom_event
phynt_set_state(const or_t3d_pos *pos, const or_t3d_att *att,
                const or_t3d_vel *vel, const or_t3d_avel *avel,
                const or_t3d_acc *acc, const or_t3d_aacc *aacc,
                const or_t3d_jerk *jerk, const or_t3d_snap *snap,
                or_rigid_body_state *reference,
                const genom_context self)
{
  (void)self; /* -Wunused-parameter */
  struct timeval tv;

  gettimeofday(&tv, NULL);
  reference->ts.sec = tv.tv_sec;
  reference->ts.nsec = tv.tv_usec * 1000.;
  reference->intrinsic = false;

  reference->pos._present = isnan(pos->x) ? false : true;
  reference->pos._value = *pos;
  reference->att._present = isnan(att->qw) ? false : true;
  reference->att._value = *att;

  reference->vel._present = isnan(vel->vx) ? false : true;
  reference->vel._value = *vel;
  reference->avel._present = isnan(avel->wx) ? false : true;
  reference->avel._value = *avel;

  reference->acc._present = isnan(acc->ax) ? false : true;
  reference->acc._value = *acc;
  reference->aacc._present = isnan(aacc->awx) ? false : true;
  reference->aacc._value = *aacc;

  reference->jerk._present = isnan(jerk->jx) ? false : true;
  reference->jerk._value = *jerk;

  reference->snap._present = isnan(snap->sx) ? false : true;
  reference->snap._value = *snap;

  return genom_ok;
}


/* --- Function set_position -------------------------------------------- */

/** Codel phynt_set_position of function set_position.
 *
 * Returns genom_ok.
 */
genom_event
phynt_set_position(double x, double y, double z, double yaw,
                   or_rigid_body_state *reference,
                   const genom_context self)
{
  (void)self; /* -Wunused-parameter */

  struct timeval tv;

  gettimeofday(&tv, NULL);
  reference->ts.sec = tv.tv_sec;
  reference->ts.nsec = tv.tv_usec * 1000.;
  reference->intrinsic = false;

  reference->pos._present = true;
  reference->pos._value.x = x;
  reference->pos._value.y = y;
  reference->pos._value.z = z;

  reference->att._present = true;
  reference->att._value.qw = cos(yaw/2.);
  reference->att._value.qx = 0.;
  reference->att._value.qy = 0.;
  reference->att._value.qz = sin(yaw/2.);

  reference->vel._present = true;
  reference->vel._value.vx = 0.;
  reference->vel._value.vy = 0.;
  reference->vel._value.vz = 0.;

  reference->avel._present = true;
  reference->avel._value.wx = 0.;
  reference->avel._value.wy = 0.;
  reference->avel._value.wz = 0.;

  reference->acc._present = true;
  reference->acc._value.ax = 0.;
  reference->acc._value.ay = 0.;
  reference->acc._value.az = 0.;

  reference->aacc._present = true;
  reference->aacc._value.awx = 0.;
  reference->aacc._value.awy = 0.;
  reference->aacc._value.awz = 0.;

  reference->jerk._present = true;
  reference->jerk._value.jx = 0.;
  reference->jerk._value.jy = 0.;
  reference->jerk._value.jz = 0.;

  reference->snap._present = true;
  reference->snap._value.sx = 0.;
  reference->snap._value.sy = 0.;
  reference->snap._value.sz = 0.;

  return genom_ok;
}


/* --- Function stop ---------------------------------------------------- */

/** Codel phynt_servo_stop of function stop.
 *
 * Returns genom_ok.
 */
genom_event
phynt_servo_stop(or_rigid_body_state *reference,
                 const genom_context self)
{
  (void)self; /* -Wunused-parameter */

  struct timeval tv;

  gettimeofday(&tv, NULL);
  reference->ts.sec = tv.tv_sec;
  reference->ts.nsec = tv.tv_usec * 1000.;
  reference->intrinsic = false;

  reference->pos._present = false;
  reference->att._present = false;
  reference->vel._present = false;
  reference->avel._present = false;
  reference->acc._present = false;
  reference->aacc._present = false;
  reference->jerk._present = false;
  reference->snap._present = false;

  return genom_ok;
}


/* --- Function log ----------------------------------------------------- */

/** Codel phynt_log of function log.
 *
 * Returns genom_ok.
 * Throws phynt_e_sys.
 */
genom_event
phynt_log(const char path[64], uint32_t decimation, phynt_log_s **log,
          const genom_context self)
{
  int fd;

  fd = open(path, O_WRONLY|O_APPEND|O_CREAT|O_TRUNC, 0666);
  if (fd < 0) return phynt_e_sys_error(path, self);

  if (write(fd, phynt_log_header_fmt "\n", sizeof(phynt_log_header_fmt)) < 0)
    return phynt_e_sys_error(path, self);

  if ((*log)->req.aio_fildes >= 0) {
    close((*log)->req.aio_fildes);

    if ((*log)->pending)
      while (aio_error(&(*log)->req) == EINPROGRESS)
        /* empty body */;
  }
  (*log)->req.aio_fildes = fd;
  (*log)->pending = false;
  (*log)->skipped = false;
  (*log)->decimation = decimation < 1 ? 1 : decimation;
  (*log)->missed = 0;
  (*log)->total = 0;

  return genom_ok;
}


/* --- Function log_stop ------------------------------------------------ */

/** Codel phynt_log_stop of function log_stop.
 *
 * Returns genom_ok.
 */
genom_event
phynt_log_stop(phynt_log_s **log, const genom_context self)
{
  (void)self; /* -Wunused-parameter */

  if (*log && (*log)->req.aio_fildes >= 0)
    close((*log)->req.aio_fildes);
  (*log)->req.aio_fildes = -1;

  return genom_ok;
}


/* --- Function log_info ------------------------------------------------ */

/** Codel phynt_log_info of function log_info.
 *
 * Returns genom_ok.
 */
genom_event
phynt_log_info(const phynt_log_s *log, uint32_t *miss, uint32_t *total,
               const genom_context self)
{
  *miss = *total = 0;
  if (log) {
    *miss = log->missed;
    *total = log->total;
  }
  return genom_ok;
}
