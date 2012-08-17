/*
 *   Copyright (c) 1996-2001 Lucent Technologies.
 *   See README file for details.
 */

#include <Python.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static PyObject *
spam_system(PyObject *self, PyObject *args)
{
    char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return Py_BuildValue("i", sts);
}
int lf_maxit = 20;
int lf_debug = 0;

static double s0, s1, tol;
static lfdata *lf_lfd;
static design *lf_des;
static smpar   *lf_sp;
int lf_status;
int ident=0;
int (*like)();
extern double robscale;

void lfdata_init(lfd)
lfdata *lfd;
{ int i;
  for (i=0; i<MXDIM; i++)
  { lfd->sty[i] = 0;
    lfd->sca[i] = 1.0;
    lfd->xl[i] = lfd->xl[i+MXDIM] = 0.0;
  }
  lfd->y = lfd->w = lfd->c = lfd->b = NULL;
  lfd->d = lfd->n = 0;
}

void smpar_init(sp,lfd)
smpar *sp;
lfdata *lfd;
{ nn(sp)  = 0.7;
  fixh(sp)= 0.0;
  pen(sp) = 0.0;
  acri(sp)= ANONE;
  deg(sp) = deg0(sp) = 2;
  ubas(sp) = 0;
  kt(sp) = KSPH;
  ker(sp) = WTCUB;
  fam(sp) = 64+TGAUS;
  link(sp)= LDEFAU;
  npar(sp) = calcp(sp,lfd->d);
}

void deriv_init(dv)
deriv *dv;
{ dv->nd = 0;
}


int locfit(lfd,des,sp,noit,nb,cv)
lfdata *lfd;
design *des;
smpar *sp;
int noit, nb, cv;
{ int i;

  if (des->xev==NULL)
  { ERROR(("locfit: NULL evaluation point?"));
    return(246);
  }

  if (lf_debug>0)
  { printf("locfit: ");
    for (i=0; i<lfd->d; i++) printf(" %10.6f",des->xev[i]);
    printf("\n");
  }

  lf_des = des;
  lf_lfd = lfd;
  lf_sp  = sp;

/* the 1e-12 avoids problems that can occur with roundoff */
  if (nb) nbhd(lfd,des,(int)(lfd->n*nn(sp)+1e-12),0,sp);

  lf_status = lfinit(lfd,sp,des);
  if (lf_status != LF_OK) return(lf_status);

  if (use_robust_scale(fam(sp)))
    lf_robust(lfd,sp,des,lf_maxit);
  else
  { robscale = 1.0;
    lfiter(des,lf_maxit);
  }

  if (lf_status == LF_OOB) setzero(des->cf,des->p);

  if ((fam(sp)&63)==TDEN) /* convert from rate to density */
  { switch(link(sp))
    { case LLOG:
        des->cf[0] -= log(des->smwt);
        break;
      case LIDENT:
        multmatscal(des->cf,1.0/des->smwt,des->p);
        break;
      default: ERROR(("Density adjustment; invalid link"));
    }
  }

  /* variance calculations, if requested */
  if (cv)
    lf_vcov(lfd,sp,des);

  return(lf_status);
}
static PyMethodDef SpamMethods[] = {
    {"spam_system", spam_system, METH_VARARGS, "Calculate the Fibonacci numbers."},
    {"locfit", locfit, METH_VARARGS, "Calculate the Fibonacci numbers."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initspam(void)
{
    (void) Py_InitModule("spam", SpamMethods);
}

