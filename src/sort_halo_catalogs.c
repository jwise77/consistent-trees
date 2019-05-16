#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>
#include <string.h>
#include "check_syscalls.h"
#include "stringparse.h"
#include "inthash.h"

#define MASS_LABEL "Mvir"

char *halo_data = NULL;
int64_t data_alloced = 0;
int64_t total_size = 0;

struct halo_stats {
  int64_t id, upid;
  float mass, host_mass;
  int64_t data_offset;
  int64_t data_length;
};

struct halo_stats *hs = NULL;
int64_t num_halos = 0;

int sort_halos(const void *a, const void *b) {
  const struct halo_stats *c = a;
  const struct halo_stats *d = b;
  if (c->host_mass > d->host_mass) return -1;
  if (d->host_mass > c->host_mass) return 1;
  if (c->upid < d->upid) return -1;
  if (d->upid < c->upid) return 1;
  if (c->mass > d->mass) return -1;
  if (d->mass > c->mass) return 1;
  return 0;
}

int main(int argc, char **argv) {
  int64_t i;
  char buffer[2048];
  if (argc < 2) {
    printf("Usage: %s /path/to/hlist\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  FILE *in = check_fopen(argv[1], "r");
  snprintf(buffer, 2048, "%s.sorted", argv[1]);
  FILE *out = check_fopen(buffer, "w");
  struct inthash *upids = new_inthash();
  struct inthash *host_masses = new_inthash();

  SHORT_PARSETYPE;
  struct halo_stats ts;
  struct parse_format pf[3] = {{1,D64,&ts.id},
			       {6,D64,&ts.upid},
			       {10,F,&ts.mass}};
  
  int64_t mass_column = 10;
  int64_t line_count = 0;
  while (fgets(buffer, 2048, in)) {
    if (buffer[0]=='#') {
      fprintf(out, "%s", buffer);
      if (line_count==0) {
	char *mass_label = " "MASS_LABEL"(";
	char *mass_str = strstr(buffer, mass_label);
	if (mass_str && sscanf(mass_str+strlen(mass_label), "%"SCNd64, &mass_column)==1) {
	  printf("Sorting by %s in column %"PRId64"\n", MASS_LABEL, mass_column);
	  pf[2].index = mass_column;
	}
      }
      line_count++;
      continue;
    }
    if (stringparse_format(buffer, pf, 3)<3) continue;
    ts.data_length = strlen(buffer);
    ts.data_offset = total_size;
    check_realloc_smart(halo_data, 1, data_alloced, (total_size+ts.data_length));
    memcpy(halo_data+ts.data_offset, buffer, ts.data_length);
    total_size += ts.data_length;
    if (ts.upid < 0) ts.upid = ts.id;
    ih_setint64(upids, ts.id, ts.upid);
    check_realloc_every(hs, sizeof(struct halo_stats), num_halos, 1000);
    hs[num_halos] = ts;
    num_halos++;
  }
  fclose(in);
    
  for (i=0; i<num_halos; i++) {
    while (1) {
      int64_t new_upid = ih_getint64(upids, hs[i].upid);
      assert(new_upid != IH_INVALID);
      if (new_upid == hs[i].upid) break;
      hs[i].upid = new_upid;
    }
    ih_setval(host_masses, hs[i].id, &(hs[i].mass));
  }

  for (i=0; i<num_halos; i++) {
    float *hm = ih_getval(host_masses, hs[i].upid);
    assert(hm != NULL);
    hs[i].host_mass = *hm;
  }
  
  qsort(hs, num_halos, sizeof(struct halo_stats), sort_halos);

  for (i=0; i<num_halos; i++)
    fwrite(halo_data+hs[i].data_offset, 1, hs[i].data_length, out);
  fclose(out);
  return 0;
}


