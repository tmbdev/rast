/* Copyright (c) 1990-1995 by Thomas M. Breuel */

/*
  Command line driver program for the RAST library.
  See the documentation for arguments and parameters.
*/

#include <stdio.h>

#include "misc.h"
#include "narray.h"
#include "smartptr.h"
#include "vec2.h"
using namespace colib;

#include "util.h"
#include "rast.h"

const char *usage =
    "Usage: rast subprogram ...\n"
    "\n"
    "Parameter settings are passed in the environment (e.g.,\n"
    "verbose=1 epsilon=5 ./rast rast model image)\n"
    "\n"
    "Use verbose_params=1 to see the settings being used.\n"
    "\n"
    "See the documentation for an explanation of what these parameters mean.\n"
    "\n"
    "rast instance model.out image.out\n"
    "   image_size model_size nclutter nmodel_unoccluded error aerror\n"
    "   minscale maxscale\n"
    "rast align model.points image.points\n"
    "   minscale maxscale error\n"
    "rast rast model.points image.points\n"
    "   maxresults verbose tolerance min_q mindx maxdx mindy maxdy\n"
    "   amin amax minscale maxscale lsq eps aeps\n"
    "rast srast model.segments image.segments\n"
    "   maxresults verbose mindx maxdx mindy maxdy amin amax minscale "
    "maxscale\n"
    "   eps aeps eps_scales ieps lsq tolerance qtolerance\n"
    "rast ssrast model.segments image.segments\n"
    "   maxresults verbose mindx maxdx mindy maxdy amin amax minscale "
    "maxscale\n"
    "   eps aeps sdist lsq tolerance qtolerance\n"
    "rast lines data.points\n"
    "   usegrad maxresults verbose error angle_error tolerance "
    "angle_tolerance\n"
    "   minweight lsq\n"
    "rast slines data.segments\n"
    "   usegrad maxresults verbose error angle_error tolerance "
    "angle_tolerance\n"
    "   minweight maxoffset lsq unoriented\n";

int main_instance(int argc, char **argv) {
  if (argc != 3) throw "wrong # args";
  srand48(igetenv("seed", mkseed()));
  autodel<InstanceP2D> instance(makeInstanceP2D());
  instance->set_image_size(igetenv("image_size", 512));
  instance->set_model_size(igetenv("model_size", 100));
  instance->set_nclutter(igetenv("nclutter", 50));
  instance->set_nmodel_total(igetenv("nmodel_total", 20));
  instance->set_nmodel_unoccluded(igetenv("nmodel_unoccluded", 10));
  float error = fgetenv("error", 5);
  float aerror = fgetenv("aerror", 0.1);
  instance->set_error(error);
  instance->set_aerror(aerror);
  instance->set_srange(fgetenv("minscale", 0.8), fgetenv("maxscale", 1.2));
  instance->generate();

  stdio model(argv[1], "w");

  for (int i = 0; i < instance->nmodel(); i++) {
    float x, y, a;
    instance->get_model(x, y, a, i);
    fprintf(model, "%g %g %g %g %g\n", x, y, a, error, aerror);
  }

  model.close();

  stdio image(argv[2], "w");

  fprintf(image, "# ");
  for (int i = 0; i < 4; i++) fprintf(image, " %g", instance->get_param(i));
  fprintf(image, "\n");

  for (int i = 0; i < instance->nimage(); i++) {
    float x, y, a;
    instance->get_image(x, y, a, i);
    fprintf(image, "%g %g %g\n", x, y, a);
  }

  image.close();
  return 0;
}

int main_align(int argc, char **argv) {
  if (argc != 3) throw "wrong # args";
  autodel<AlignmentP2D> align(makeAlignmentP2D());
  align->set_srange(fgetenv("minscale", 0.8), fgetenv("maxscale", 1.2));
  align->set_epsilon(fgetenv("error", 5.0));
  char buf[1000];
  stdio model(argv[1], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, model)) break;
    if (buf[0] == '#') continue;
    float x, y;
    if (sscanf(buf, "%g %g", &x, &y) != 2) throw "bad format";
    align->add_mpoint(x, y);
  }
  stdio image(argv[2], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, image)) break;
    if (buf[0] == '#') continue;
    float x, y;
    if (sscanf(buf, "%g %g", &x, &y) != 2) throw "bad format";
    align->add_ipoint(x, y);
  }
  align->compute();
  printf("%g    %g %g %g %g\n", align->quality(), align->translation(0),
         align->translation(1), align->angle(), align->scale());
  return 0;
}

int main_rastp2d(int argc, char **argv) {
  if (argc != 3) throw "wrong # args";
  autodel<RastP2D> rast(makeRastP2D());
  rast->set_maxresults(igetenv("maxresults", 1));
  rast->set_verbose(igetenv("verbose", 0));
  rast->set_tolerance(fgetenv("tolerance", 1e-3));
  rast->set_min_q(fgetenv("min_q", 3));
  rast->set_xrange(fgetenv("mindx", -1000), fgetenv("maxdx", 1000));
  rast->set_yrange(fgetenv("mindy", -1000), fgetenv("maxdy", 1000));
  rast->set_arange(fgetenv("amin", 0.0), fgetenv("amax", 2 * M_PI));
  rast->set_srange(fgetenv("minscale", 0.8), fgetenv("maxscale", 1.2));
  rast->set_lsq(igetenv("lsq", 0));
  float eps = fgetenv("eps", 5.0);
  float aeps = fgetenv("aeps", 0.1);
  char buf[1000];
  stdio model(argv[1], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, model)) break;
    if (buf[0] == '#') continue;
    float x, y, a, err = eps, aerr = aeps;
    int nfields = sscanf(buf, "%g %g %g %g %g", &x, &y, &a, &err, &aerr);
    if (nfields < 3) throw "bad format";
    rast->add_msource(x, y, a, err, aerr);
  }
  stdio image(argv[2], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, image)) break;
    if (buf[0] == '#') continue;
    float x, y, a;
    if (sscanf(buf, "%g %g %g", &x, &y, &a) != 3) throw "bad format";
    rast->add_ipoint(x, y, a);
  }
  rast->match();
  for (int i = 0; i < rast->nresults(); i++) {
    printf("%d  %g %g   %g %g %g %g\n", i, rast->ubound(i), rast->lbound(i),
           rast->translation(i, 0), rast->translation(i, 1), rast->angle(i),
           rast->scale(i));
  }
  rast.move();  // don't waste time deallocating
  return 0;
}

int main_rasts2d(int argc, char **argv) {
  if (argc != 3) throw "wrong # args";
  autodel<RastS2D> rast(makeRastS2D());
  rast->set_maxresults(igetenv("maxresults", 1));
  rast->set_verbose(igetenv("verbose", 0));
  rast->set_xrange(fgetenv("mindx", -1000), fgetenv("maxdx", 1000));
  rast->set_yrange(fgetenv("mindy", -1000), fgetenv("maxdy", 1000));
  rast->set_arange(fgetenv("amin", 0.0), fgetenv("amax", 2 * M_PI));
  rast->set_srange(fgetenv("minscale", 0.8), fgetenv("maxscale", 1.2));
  float eps = fgetenv("eps", 4.0);
  float aeps = fgetenv("aeps", 0.05);
  rast->set_eps(eps, aeps);
  rast->set_scale_eps(igetenv("eps_scales", 0), fgetenv("ieps", eps));
  rast->set_lsq(igetenv("lsq", 0));
  rast->set_tolerance(fgetenv("tolerance", 1e-3));
  rast->set_qtolerance(fgetenv("qtolerance", 1e-2));
  //!!!float min_q_frac = fgetenv("min_q_frac",0.5);
  //!!!rast->set_match_mode(sgetenv("match_mode","sub_lsq"));
  char buf[1000];
  stdio model(argv[1], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, model)) break;
    if (buf[0] == '#') continue;
    float x, y, x1, y1;
    if (sscanf(buf, "%g %g %g %g", &x, &y, &x1, &y1) < 4) throw "bad format";
    rast->add_mseg(x, y, x1, y1);
  }
  stdio image(argv[2], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, image)) break;
    if (buf[0] == '#') continue;
    float x, y, x1, y1;
    if (sscanf(buf, "%g %g %g %g", &x, &y, &x1, &y1) < 4) throw "bad format";
    rast->add_iseg(x, y, x1, y1);
  }
  // rast->set_min_q(rast->total_weight()*min_q_frac);
  double start = now();
  rast->match();
  double total = now() - start;
  fprintf(stderr, "time %g\n", total);
  for (int i = 0; i < rast->nresults(); i++) {
    printf("%d  %g %g   %g %g %g %g\n", i, rast->ubound(i), rast->lbound(i),
           rast->translation(i, 0), rast->translation(i, 1), rast->angle(i),
           rast->scale(i));
  }
  fflush(stdout);
  rast.move();
  return 0;
}

int main_rastss2d(int argc, char **argv) {
  if (argc != 3) throw "wrong # args";
  autodel<RastSS2D> rast(makeRastSS2D());
  rast->set_maxresults(igetenv("maxresults", 1));
  rast->set_verbose(igetenv("verbose", 0));
  rast->set_xrange(fgetenv("mindx", -1000), fgetenv("maxdx", 1000));
  rast->set_yrange(fgetenv("mindy", -1000), fgetenv("maxdy", 1000));
  rast->set_arange(fgetenv("amin", 0.0), fgetenv("amax", 2 * M_PI));
  rast->set_srange(fgetenv("minscale", 0.8), fgetenv("maxscale", 1.2));
  float eps = fgetenv("eps", 4.0);
  float aeps = fgetenv("aeps", 0.05);
  float sdist = fgetenv("sdist", eps);
  rast->set_eps(eps, aeps, sdist);
  rast->set_lsq(igetenv("lsq", 0));
  rast->set_tolerance(fgetenv("tolerance", 1e-3));
  rast->set_qtolerance(fgetenv("qtolerance", 1e-2));
  //!!!float min_q_frac = fgetenv("min_q_frac",0.5);
  //!!!rast->set_match_mode(sgetenv("match_mode","sub_lsq"));
  char buf[1000];
  stdio model(argv[1], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, model)) break;
    if (buf[0] == '#' || buf[0] == '\0' || buf[0] == '\n') continue;
    float x, y, x1, y1;
    if (sscanf(buf, "%g %g %g %g", &x, &y, &x1, &y1) < 4) continue;
    rast->add_mseg(x, y, x1, y1);
  }
  stdio image(argv[2], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, image)) break;
    if (buf[0] == '#' || buf[0] == '\0' || buf[0] == '\n') continue;
    float x, y, x1, y1;
    if (sscanf(buf, "%g %g %g %g", &x, &y, &x1, &y1) < 4) throw "bad format";
    rast->add_iseg(x, y, x1, y1);
  }
  // rast->set_min_q(rast->total_weight()*min_q_frac);
  double start = now();
  rast->match();
  double total = now() - start;
  fprintf(stderr, "time %g\n", total);
  for (int i = 0; i < rast->nresults(); i++) {
    printf("%d  %g %g   %g %g %g %g\n", i, rast->ubound(i), rast->lbound(i),
           rast->translation(i, 0), rast->translation(i, 1), rast->angle(i),
           rast->scale(i));
  }
  fflush(stdout);
  rast.move();
  return 0;
}

int main_lines(int argc, char **argv) {
  if (argc != 2) throw "wrong # args";
  autodel<LinesP2D> rast(makeLinesP2D());
  bool usegrad = igetenv("usegrad", 0);
  bool useweights = igetenv("useweights", 0);
  rast->set_maxresults(igetenv("maxresults", 1));
  rast->set_verbose(igetenv("verbose", 0));
  double angle_error = fgetenv("angle_error", 0.1);
  rast->set_error(fgetenv("error", 2.0), angle_error);
  rast->set_tolerance(fgetenv("tolerance", 0.1),
                      fgetenv("angle_tolerance", 0.001));
  rast->set_minweight(fgetenv("minweight", 0.0));
  rast->set_lsq(igetenv("lsq", 1));
  char buf[1000];
  stdio points(argv[1], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, points)) break;
    if (buf[0] == '#') continue;
    if (buf[0] == '\n') continue;
    float x, y, a = 0, w = 1;
    int fields = sscanf(buf, "%g %g %g %g", &x, &y, &a, &w);
    if (usegrad) a = normangleOf(a - M_PI / 2.0);
    if (fields < 2 && angle_error < 2 * M_PI)
      throw "no angles available; set angle_error to >=2pi";
    if (fields < 2 || fields > 4) throw "bad format";
    if (!useweights) w = 1.0;
    rast->add_ipoint(x, y, a, w);
  }
  rast->compute();
  for (int i = 0; i < rast->nresults(); i++) {
#if 1
    printf("%d  %g   %g %g\n", i, rast->weight(i), rast->angle(i),
           rast->offset(i));
#else
    printf("rank %d quality %g nmatches %d params %g %g\n", i, rast->weight(i),
           rast->nmatches(i), rast->angle(i), rast->offset(i));
#endif
  }
  return 0;
}

int main_slines(int argc, char **argv) {
  if (argc != 2) throw "wrong # args";
  autodel<LinesS2D> rast(makeLinesS2D());
  bool usegrad = igetenv("usegrad", 0);
  rast->set_maxresults(igetenv("maxresults", 1));
  rast->set_verbose(igetenv("verbose", 0));
  double angle_error = fgetenv("angle_error", 0.1);
  rast->set_error(fgetenv("error", 2.0), angle_error);
  rast->set_tolerance(fgetenv("tolerance", 0.1),
                      fgetenv("angle_tolerance", 0.001));
  rast->set_minweight(fgetenv("minweight", 0.0));
  rast->set_maxoffset(fgetenv("maxoffset", 3000.0));
  rast->set_lsq(igetenv("lsq", 1));
  rast->set_unoriented(igetenv("unoriented", 1));
  char buf[1000];
  stdio points(argv[1], "r");
  for (;;) {
    if (!fgets(buf, sizeof buf, points)) break;
    if (buf[0] == '#') continue;
    if (buf[0] == '\n') continue;
    float x, y, x1, y1, a = 0, w = 0;
    int fields = sscanf(buf, "%g %g %g %g %g %g", &x, &y, &x1, &y1, &a, &w);
    if (usegrad) a = normangleOf(a - M_PI / 2.0);
    if (fields < 4 && angle_error < 2 * M_PI)
      throw "no angles available; set angle_error to >=2pi";
    if (fields < 4 || fields > 6) throw "bad format";
    if (w == 0) w = sqrt(sqr(x1 - x) + sqr(y1 - y));
    rast->add_iseg(x, y, x1, y1, a, w);
  }
  double start = now();
  rast->compute();
  double total = now() - start;
  fprintf(stderr, "time %g\n", total);
  for (int i = 0; i < rast->nresults(); i++) {
    printf("%d  %g   %g %g\n", i, rast->weight(i), rast->angle(i),
           rast->offset(i));
  }
  fprintf(stderr, "# done\n");
  rast.move();
  return 0;
}

int main(int argc, char **argv) {
  try {
    if (argc < 2) {
      fprintf(stderr, "%s", usage);
      exit(255);
    }
    if (!strcmp(argv[1], "instance")) return main_instance(argc - 1, argv + 1);
    if (!strcmp(argv[1], "align")) return main_align(argc - 1, argv + 1);
    if (!strcmp(argv[1], "lines")) return main_lines(argc - 1, argv + 1);
    if (!strcmp(argv[1], "slines")) return main_slines(argc - 1, argv + 1);
    if (!strcmp(argv[1], "rast")) return main_rastp2d(argc - 1, argv + 1);
    if (!strcmp(argv[1], "srast")) return main_rasts2d(argc - 1, argv + 1);
    if (!strcmp(argv[1], "ssrast")) return main_rastss2d(argc - 1, argv + 1);
    fprintf(stderr, "unknown subprogram\n");
    exit(1);
  } catch (const char *error) {
    fprintf(stderr, "error: %s\n", error);
    exit(2);
  }
}
