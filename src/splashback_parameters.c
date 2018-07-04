#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "halo_io.h"
#include "grav_config.h"
#include "gravitational_consistency_vars.h"
#include "check_syscalls.h"
#include "masses.h"
#include "halo_io.h"
#include "inthash.h"
#include "version.h"
#include "stringparse.h"

struct halo_stash now={0}, evolved = {0}, dead = {0};
float box_size=0;
float max_mvir=0;
float min_mvir=0;
int64_t *sort_order = NULL;
struct inthash *forest=NULL, *forest_evolved=NULL, *forest_now=NULL;

#define MAX_MASSES 10
int64_t num_masses = 0;
char *mass_defs[MAX_MASSES] = {0};
double r_convs[MAX_MASSES] = {0};
char *header_line = NULL;
int64_t mass_columns[MAX_MASSES] = {0};

struct infall_stats {
  int64_t id, descid, upid, infall_desc_upid;
  int64_t flags;
  int32_t mmp;
  float pos[6], m, r, vmax;
  float infall_scale, infall_mass, infall_vmax, infall_peak_mass, infall_rfrac, peak_mass, peak_vmax, v_tan, v_rad, vhost, mhost, rhost;
};

struct infall_stats *nis = NULL;
struct infall_stats *ois = NULL;
struct inthash *ni = NULL, *oi = NULL;

int64_t num_new = 0;
int64_t num_old = 0;


int strcmp_sort(const void *a, const void *b) {
  char * const *c = a;
  char * const *d = b;
  assert(strncmp(*c, "hlist_", 6)==0);
  assert(strncmp(*d, "hlist_", 6)==0);
  double e = atof(*c+6);
  double f = atof(*d+6);
  if (e<f) return -1;
  if (f<e) return 1;
  return 0;
}

void populate_host_props(struct infall_stats *infall, struct infall_stats *target, struct infall_stats *host, float scale) {
  int64_t i;
  float pos[6], vt=0, vr=0, ds=0;
  for (i=0; i<6; i++) {
    pos[i] = infall->pos[i] - host->pos[i];
    if (i<3) {
      if (pos[i] > box_size/2.0) pos[i] -= box_size;
      else if (pos[i] < -box_size/2.0) pos[i] += box_size;
    }
  }
  for (i=0; i<3; i++) ds += pos[i]*pos[i];
  ds = sqrt(ds);
  if (ds>0) for (i=0; i<3; i++) pos[i]/=ds;
  for (i=0; i<3; i++) {
    float dvr = pos[i+3]*pos[i];
    vr += dvr;
  }
  for (i=0; i<3; i++) {
    float dvt = pos[i+3] - vr*pos[i];
    vt += dvt*dvt;
  }
  vt = sqrt(vt);
  target->infall_rfrac = ds * 1e3 / host->r;
  target->v_tan = vt;
  target->v_rad = vr;
  target->mhost = host->m;
  target->rhost = host->r;
  target->vhost = 2.07419097e-3*sqrt(host->m / (scale*host->r)); //sqrt(G*solar mass / kpc) = 2.07... m/s
}

int main(int argc, char **argv) {
  int64_t i, j;
  char buffer[2048];
  float scale = 0, prev_scale = 0;
  
  if (argc==1) {
    fprintf(stderr, "Consistent Trees, Version %s\n", TREE_VERSION);
    fprintf(stderr, "%s.  See the LICENSE file for redistribution details.\n", TREE_COPYRIGHT);
    fprintf(stderr, "Usage: %s hlist*.list \n", argv[0]); exit(EXIT_FAILURE);
  }

  qsort(argv+1, argc-1, sizeof(char *), strcmp_sort);

  oi = new_inthash();
  
  for (i=1; i<argc; i++) {
    for (j=0; argv[i][j]!=0 && (argv[i][j]<'0' || argv[i][j]>'9'); j++);
    if (argv[i][j]==0) {
      fprintf(stderr, "Failed to parse scale factor in filename \"%s\"\n", argv[i]);
      exit(EXIT_FAILURE);
    }
    scale = atof(argv[i]+j);

    FILE *in = check_fopen(argv[i], "r");
    snprintf(buffer, 2048, "%s.splashback_info", argv[i]);
    FILE *out = check_fopen(buffer, "w");
    struct infall_stats is = {0};
    is.infall_desc_upid = -1;
    SHORT_PARSETYPE;
    #define NUM_FIELDS 13
    //scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26)
    struct parse_format pf[NUM_FIELDS] = {{1,D64,&is.id}, {3,D64,&is.descid}, {6, D64, &is.upid}, {10, F, &is.m}, {11, F, &is.r}, {14, D, &is.mmp}, {16, F, &is.vmax}, {17, F, is.pos}, {18, F, is.pos+1}, {19, F, is.pos+2}, {20, F, is.pos+3}, {21, F, is.pos+4}, {22, F, is.pos+5}};
    fgets(buffer, 2048, in);
    fprintf(out, "#ID DescID UPID Infall_UPID_Desc Flags Mass Radius Vmax X Y Z VX VY VZ Infall_Scale Infall_Mass Infall_Vmax Infall_R/R_host Peak_Mass Peak_Vmax V_Tan V_Rad Vhost Mhost Rhost\n");
    fprintf(out, "#Infall_UPID_Desc: descendant id of infall halo's upid at the current scale factor\n");
    fprintf(out, "#Flags: 1 if Infall_UPID_Desc has ever been != to current halo's UPID (-1 OK); 0 otherwise.\n");
    while (fgets(buffer, 2048, in)) {
      if (buffer[0]=='#') {
	if (!strncmp(buffer, "#Full box size = ", 17))
	  box_size = atof(buffer+17);
	fprintf(out, "%s", buffer);
	continue;
      }
      if (stringparse_format(buffer, pf, NUM_FIELDS)!=NUM_FIELDS) continue;

      is.peak_mass = is.m;
      is.peak_vmax = is.vmax;
      check_realloc_every(nis, sizeof(struct infall_stats), num_new, 1000);
      nis[num_new] = is;
      num_new++;
    }    
    fclose(in);

    struct inthash *mmp = new_inthash();
    ni = new_inthash();
    for (j=0; j<num_new; j++) ih_setval(ni, nis[j].id, nis+j);
    for (j=0; j<num_old; j++) if (ois[j].mmp) ih_setval(mmp, ois[j].descid, ois+j);
    for (j=0; j<num_new; j++) {
      struct infall_stats *prog = ih_getval(mmp, nis[j].id);

      if (!prog) {
	if (nis[j].upid > -1) { //Pristine Subhalo
	  struct infall_stats *parent = ih_getval(ni, nis[j].upid);
	  assert(parent);
	  nis[j].infall_scale = scale;
	  nis[j].infall_mass = nis[j].m;
	  nis[j].infall_peak_mass = nis[j].m;
	  nis[j].infall_vmax = nis[j].vmax;
	  nis[j].infall_desc_upid = nis[j].upid;
	  populate_host_props(nis+j, nis+j, parent, scale);
	}
	continue;
      }

      int64_t keep_infall_stats = 0;
      if (prog->infall_scale && prog->infall_peak_mass*2.0 > nis[j].m) keep_infall_stats = 1;
      
      if (prog->upid < 0 && nis[j].upid > -1) { //Just fell in
	if (!prog->infall_scale || !keep_infall_stats) {
	  nis[j].infall_scale = prev_scale;
	  nis[j].infall_mass = prog->m;
	  nis[j].infall_peak_mass = prog->peak_mass;
	  nis[j].infall_vmax = prog->vmax;
	  nis[j].infall_desc_upid = nis[j].upid;
	  struct infall_stats *parent = ih_getval(mmp, nis[j].upid);
	  if (parent) populate_host_props(prog, nis+j, parent, prev_scale);
	  else {
	    struct infall_stats *host = ih_getval(ni, nis[j].upid);
	    assert(host);
	    populate_host_props(nis+j, nis+j, host, scale);
	  }
	}
      } else if (nis[j].upid > -1) { //Remains a subhalo
	assert(prog->infall_scale);
	keep_infall_stats = 1;
      }
	
      if (prog->infall_scale && keep_infall_stats) {
	nis[j].infall_scale = prog->infall_scale;
	nis[j].infall_mass = prog->infall_mass;
	nis[j].infall_peak_mass = prog->infall_peak_mass;
	nis[j].infall_vmax = prog->infall_vmax;
	nis[j].infall_rfrac = prog->infall_rfrac;
	nis[j].v_tan = prog->v_tan;
	nis[j].v_rad = prog->v_rad;
	nis[j].vhost = prog->vhost;
	nis[j].mhost = prog->mhost;
	nis[j].rhost = prog->rhost;
	nis[j].flags = prog->flags;
	struct infall_stats *parent_prog = ih_getval(oi, prog->infall_desc_upid);
	assert(parent_prog);
	struct infall_stats *parent_now = ih_getval(ni, parent_prog->descid);
	nis[j].infall_desc_upid = parent_now->id;
	if (nis[j].upid > -1 && (parent_now->id != nis[j].upid)) nis[j].flags = 1;
      }

      if (nis[j].m < prog->peak_mass) nis[j].peak_mass = prog->peak_mass;
      if (nis[j].vmax < prog->peak_vmax) nis[j].peak_vmax = prog->peak_vmax;

      fprintf(out, "%"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId64" %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      nis[j].id, nis[j].descid, nis[j].upid, nis[j].infall_desc_upid, nis[j].flags, nis[j].m, nis[j].r, nis[j].vmax, nis[j].pos[0], nis[j].pos[1], nis[j].pos[2], nis[j].pos[3], nis[j].pos[4], nis[j].pos[5], nis[j].infall_scale, nis[j].infall_mass, nis[j].infall_vmax, nis[j].infall_rfrac, nis[j].peak_mass, nis[j].peak_vmax, nis[j].v_tan, nis[j].v_rad, nis[j].vhost, nis[j].mhost, nis[j].rhost);
    }

    free_inthash(mmp);
    if (oi) free_inthash(oi);
    oi = ni; ni = NULL;
    if (ois) free(ois);
    ois = nis;
    nis = NULL;
    num_old = num_new;
    num_new = 0;
    fclose(out);
    prev_scale = scale;
  }
  
  return 0;
}

void parse_header_line(char *hl) {
  int64_t i = 0, column=0;
  num_masses = 0;  
  header_line = strdup(hl);
  for (; header_line[i] && num_masses < MAX_MASSES; i++) {
    if (header_line[i] == '#') continue;
    if (header_line[i] == ' ') {
      if (i==0 || header_line[i-1] != ' ') column++;
      continue;
    }

    if ((header_line[i] == 'm') || (header_line[i] == 'M')) {
      i++;
      if ((header_line[i] == 'm') || (header_line[i] == 'M')) i++;
      if (strncasecmp(header_line, "vir", 3)==0) {
	if (header_line[i+3] != ' ' && header_line[i+3] != '(') continue;
	mass_defs[num_masses] = header_line+i;
	mass_columns[num_masses] = column;
	i+=3;
	header_line[i] = 0;
	num_masses++;
	continue;
      }

      if (header_line[i] >= '0' && header_line[i] <= '9') {
	int64_t j=i;
	while (header_line[j+1] != ' ' && header_line[j+1] != 0) j++;
	if (header_line[j]==')') { while (j>i && header_line[j]!= '(') j--; }
	mass_defs[num_masses] = header_line+i;
	mass_columns[num_masses] = column;
	i=j;
	header_line[i] = 0;
	num_masses++;
	continue;
      }
    }
  }
}

