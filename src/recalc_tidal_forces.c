#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <inttypes.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "gravitational_consistency.h"
#include "gravitational_consistency_vars.h"
#include "grav_config.h"
#include "distance.h"
#include "universe_time.h"
#include "check_syscalls.h"
#include "halo_io.h"
#include "version.h"
#include "tidal_lib.h"

int64_t MAX_FORKS=6;

struct halo_stash now={0};
int64_t children = 0;

float box_size = 0; /* Automatically set; in comoving Mpc/h */
float max_mvir = 0; /* Automatically set; in Msun */
float min_mvir = 0;

void clear_halo_stash(struct halo_stash *h) {
  if (h->halos) h->halos = (struct tree_halo *)realloc(h->halos, 0);
  if (h->id_conv) h->id_conv = (int64_t *)realloc(h->id_conv, 0);
  h->max_id = h->min_id = h->num_halos = 0;
  //max_mvir = 0;
}

void wait_for_children(int all) {
  int stat_loc;
  for (; children >= MAX_FORKS; children--)
    if (!wait4(-1, &stat_loc, 0, 0))
      children = 1;
  if (all) {
    while (wait4(-1, &stat_loc, 0, 0)>0);
    children = 0;
  }
}

int main(int argc, char **argv) {
  int64_t i, num_outputs=0;
  float *output_scales=NULL;
  int64_t *outputs=NULL;
  pid_t pid;
  float a1;
  char buffer[1024];

  if (argc==1) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s options.cfg\n", argv[0]); exit(1);
  }
  grav_config(argv[1], 0);
  print_timing(_RTF, NULL);

  init_cosmology(Om, Ol, h0);
  init_time_table(Om, h0);
  read_outputs(&output_scales, &outputs, &num_outputs);

  for (i=0; i<num_outputs; i++)
  {
    a1 = output_scales[i];
    timed_output(_RTF, "** Starting work for snapshot %"PRId64" (a=%f)...", outputs[i], a1);

    clear_halo_stash(&now);
    timed_output(_RTF, "Loading snapshot %"PRId64"...", outputs[i]);
    snprintf(buffer, 1024, "%s/consistent_%"PRId64".list", OUTBASE, outputs[i]);
    load_halos(buffer, &now, a1, 0);
    build_id_conv_list(&now);

    tidal_extra_range = tidal_fpac_mode = 1;
    timed_output(_RTF, "Calculating tidal forces...");
    calc_tidal_forces(&now, a1, a1);
    timed_output(_RTF, "Rescaling tidal forces...");
    scale_tidal_forces(&now, a1, a1);

    timed_output(_RTF, "Forking and writing halos for snapshot %"PRId64"...", outputs[i]);
    pid = fork();
    if (pid < 1) {
      char new_filename[1024];
      snprintf(buffer, 1024, "%s/really_consistent_%"PRId64".list", OUTBASE, outputs[i]);
      snprintf(new_filename, 1024, "%s/really_consistent_%"PRId64".list.new", OUTBASE, outputs[i]);
      common_print_halos(new_filename, now.halos, now.num_halos, 0);
      rename(new_filename, buffer);
      if (!pid) return 0;
    } else {
      children++;
    }
    wait_for_children(0);
    timed_output(_RTF, "** Done with snapshot %"PRId64".", outputs[i]);
  }
  wait_for_children(1);
  timed_output(_RTF, "Successfully finished.");
  close_timing_log();
  return 0;
}


