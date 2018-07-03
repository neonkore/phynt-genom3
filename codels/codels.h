/*
 * Copyright (c) 2016-2018 LAAS/CNRS
 * All rights reserved.
 *
 * Redistribution  and  use  in  source  and binary  forms,  with  or  without
 * modification, are permitted provided that the following conditions are met:
 *
 *   1. Redistributions of  source  code must retain the  above copyright
 *      notice and this list of conditions.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice and  this list of  conditions in the  documentation and/or
 *      other materials provided with the distribution.
 *
 * THE SOFTWARE  IS PROVIDED "AS IS"  AND THE AUTHOR  DISCLAIMS ALL WARRANTIES
 * WITH  REGARD   TO  THIS  SOFTWARE  INCLUDING  ALL   IMPLIED  WARRANTIES  OF
 * MERCHANTABILITY AND  FITNESS.  IN NO EVENT  SHALL THE AUTHOR  BE LIABLE FOR
 * ANY  SPECIAL, DIRECT,  INDIRECT, OR  CONSEQUENTIAL DAMAGES  OR  ANY DAMAGES
 * WHATSOEVER  RESULTING FROM  LOSS OF  USE, DATA  OR PROFITS,  WHETHER  IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR  OTHER TORTIOUS ACTION, ARISING OUT OF OR
 * IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *                                           Anthony Mallet on Tue Jun 12 2018
 */
#ifndef H_PHYNT_CODELS
#define H_PHYNT_CODELS

#include <aio.h>
#include <errno.h>
#include <string.h>

#include "phynt_c_types.h"

#ifdef __cplusplus
extern "C" {
#endif

  int	phynt_admittance_filter(const phynt_ids_body_s *body,
                                const phynt_ids_af_s *af,
                                const or_rigid_body_state *reference,
                                const or_wrench_estimator_state *ewrench,
                                or_rigid_body_state *desired,
                                phynt_log_s *log);
  void	phynt_admittance_J(const double J[3 * 3]);

  int	phynt_wrench_observer(const phynt_ids_body_s *body,
                              const phynt_ids_wo_s *wo,
                              const or_pose_estimator_state *state,
                              const or_wrench_estimator_state *mwrench,
                              or_wrench_estimator_state *ewrench);

#ifdef __cplusplus
}
#endif

static inline genom_event
phynt_e_sys_error(const char *s, genom_context self)
{
  phynt_e_sys_detail d;
  size_t l = 0;

  d.code = errno;
  if (s) {
    strncpy(d.what, s, sizeof(d.what) - 3);
    l = strlen(s);
    strcpy(d.what + l, ": ");
    l += 2;
  }
  if (strerror_r(d.code, d.what + l, sizeof(d.what) - l)) {
    /* ignore error*/;
  }
  return phynt_e_sys(&d, self);
}

struct phynt_log_s {
  struct aiocb req;
  char buffer[4096];
  bool pending, skipped;
  uint32_t decimation;
  size_t missed, total;

# define phynt_g	" %g "
# define phynt_log_header_fmt                                           \
  "ts "                                                                 \
  "efx efy efz wtx ety etz "                                            \
  "x y z roll pitch yaw "                                               \
  "vx vy vz wx wy wz "                                                  \
  "ax ay az"
# define phynt_log_fmt                                                  \
  "%d.%09d "                                                            \
  phynt_g phynt_g phynt_g phynt_g phynt_g phynt_g                       \
  phynt_g phynt_g phynt_g phynt_g phynt_g phynt_g                       \
  phynt_g phynt_g phynt_g phynt_g phynt_g phynt_g                       \
  phynt_g phynt_g phynt_g
};

#endif /* H_PHYNT_CODELS */
