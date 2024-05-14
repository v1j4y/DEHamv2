#include <stdio.h>
#include <petsctime.h>
#include <slepceps.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "doubleexchange.h"
#include "get_proj_general.h"

#define W_2holes 0.8685623084154933 
#define W_3holes 0.7263427261371458
#define W_4holes 0.5962122942961313
#define W_5holes 0.4839876699122871 
#define W_6holes 

/*
 * 
 *-------------------------------------------
 * Calculate the projection for general size
 *-------------------------------------------
 * Input
 * =====
 * valxr    = The full vector
 * Istart   = Local starting id
 * Iend     = Local ending id
 * Output
 * =====
 * projvec
 */
void get_proj_general(Vec valxr, 
                      PetscInt *Istart, 
                      PetscInt *Iend, 
                      size_t *natom, 
                      int iroot, 
                      double *projvec, 
                      double *projvec2, 
                      double *normproj, 
                      const int natomax,
                      int sze,
                      int MS2,
                      int nholes,
                      double XS,
                      double V,
                      int colm,
                      //double *projmatrix,
                      size_t *configList,
                      size_t sizeCFG,
                      size_t *csfList,
                      size_t sizeCSF,
                      const igraph_t* graph){

  int			 ideter[natomax];
  int			 ideter2[natomax];
  int 		   	 kk,kko,kok,kkio;
  int 		   	 mmo,mom,mmio;
  int 		   	 p, q;
  long int       ii;
  PetscInt      iiii;
  long int       iii;
  long int       iaa2, iaa;
  long int       nrow=-1, ncol=-1;
  double rest=0.0, norm=0.0;
  int sumMs=0;
  int ntrouGauch = 0;
  int ntrouDroit = 0;
  int idx, idxprojm;
  int checkIonic;
  int pow3 = (int)pow(3,nholes);
  size_t nelec = 2*(*natom) - nholes;
  size_t nsites = natom[0];

  // Create a new dictionary
  struct Dictionary *d = dictionary_new();
  struct Dictionary *dadr = dictionary_new();
  
  idx = prepare_dictionary(d, sze, MS2);
  //printf(" Done first dict %d \n",idx);
  idx = prepare_dictionary_for_adressing(dadr, sze, MS2);
  //printf(" Done second dict %d \n",idx);
  //printf(" %s \n",dictionary_get(dadr,intToCharPtr(5)));

  (*normproj) = 0.0;

  int dim;
  dim = pow3*idx;
  double projvecmat[dim];
  double projvec3[idx];
  double projvec4[idx];
  for(ii=0;ii<dim;++ii) projvecmat[ii] = 0.0;
  for(ii=0;ii<idx;++ii) projvec[(iroot)*idx + ii] = 0.0;
  for(ii=0;ii<idx;++ii) projvec2[(iroot)*idx + ii] = 0.0;
  for(ii=0;ii<idx;++ii) projvec3[ii] = 0.0;
  for(ii=0;ii<idx;++ii) projvec4[ii] = 0.0;
  double msloc;
  //double normvec[dim];
  int ms[nholes];
  int idh[nholes];

  PetscInt ix[1];
  PetscScalar wfy[1];

  if(fabs(XS-0.0) < 1E-10){
    idxprojm = 0;
  }
  else if(fabs(XS-0.5) < 1E-10){
    idxprojm = 0;
  }
  else if(fabs(XS-1.0) < 1E-10){
    idxprojm = 1;
  }
  else if(fabs(XS-1.5) < 1E-10){
    idxprojm = 1;
  }
  else if(fabs(XS-2.0) < 1E-10){
    idxprojm = 2;
  }
  else if(fabs(XS-2.5) < 1E-10){
    idxprojm = 2;
  }
  else if(fabs(XS-3.0) < 1E-10){
    idxprojm = 3;
  }
  else if(fabs(XS-3.5) < 1E-10){
    idxprojm = 3;
  }
  else if(fabs(XS-4.0) < 1E-10){
    idxprojm = 4;
  }
  else if(fabs(XS-4.5) < 1E-10){
    idxprojm = 4;
  }
  else if(fabs(XS-5.0) < 1E-10){
    idxprojm = 5;
  }
  else if(fabs(XS-5.5) < 1E-10){
    idxprojm = 5;
  }
  else if(fabs(XS-6.0) < 1E-10){
    idxprojm = 6;
  }
  else if(fabs(XS-6.5) < 1E-10){
    idxprojm = 6;
  }
  else if(fabs(XS-7.0) < 1E-10){
    idxprojm = 7;
  }
  else if(fabs(XS-7.5) < 1E-10){
    idxprojm = 7;
  }
  else if(fabs(XS-8.0) < 1E-10){
    idxprojm = 8;
  }
  else if(fabs(XS-8.5) < 1E-10){
    idxprojm = 8;
  }
  // 9_3h ms=6.5
  if     (fabs(XS-7.5) < 1E-10){ 
    idxprojm = 2;
  }
  else if(fabs(XS-6.5) < 1E-10){
    idxprojm = 1;
  }
  else if(fabs(XS-6.5) < 1E-10){
    idxprojm = 0;
  }

  int sizefac = (int)pow(3,nholes);
  double fachole[sizefac];
  //printf(" idx=%d dim=%d sizefac=%d sze=%d MS2=%d\n",idx,dim,sizefac,sze,MS2);
  V = 0.0;
  prepareHueckelFactors(nholes, fachole, V, sizefac);

  size_t icfg[1];
  size_t icsf[1];
  size_t cfgid, csfid;

  for(ii=*Istart;ii<*Iend;ii++) {
    iii = ii + 1;
    iiii = ii;
    size_t posi = ii;
    cfgid = findCFGID(posi, sizeCFG, sizeCSF);
    csfid = findCSFID(posi, sizeCFG, sizeCSF);
    icfg[0] = configList[cfgid];
    icsf[0] = csfList[csfid];
    //printf(" %d %d | %ld %ld \n",cfgid, csfid, icfg[0], icsf[0]);
    getdetAlphaBeta (icfg[0], icsf[0], ideter, nsites, nelec);
    //printf("      \n");
    //printf("%d -- \n",ii);
    //for(int i1=0;i1<2*nsites;++i1){
    //  printf(" %d ",ideter[i1]);
    //}
    //printf("\n");
    //getdet_(&iii, ideter);
    //normproj += valxr[iiii]*valxr[iiii];
    
    checkIonic = 0;
    // Set ms to 0
    for(kk=0;kk<nholes;++kk)
      ms[kk]=0;
    
    // Set holes to -1
    for(kk=0;kk<nholes;++kk)
      idh[kk]=-1;
   
    // Calculate MS
      //printf(" Calc MS\n");
    for(kk=0;kk<nholes;++kk){
      // Find position of hole
      for(kko=kk*3+0;kko<(kk+1)*3;++kko){
        if(ideter[kko]==3){
          idh[kk] = kko - kk*3;
        }
      }
      msloc = 0;
      for(kko=kk*3+0;kko<(kk+1)*3;++kko){
        if(ideter[kko]==1){
          msloc += 1;
        }
      }
      //printf(" %5.3f ",msloc);
      for(kko=(nholes*3*2)-1-kk*3;kko>(nholes*3*2)-1-(kk+1)*3;--kko){
        if(ideter[kko]==1){
          msloc += 1;
        }
      }
      //printf(" %5.3f ",msloc);
      ms[kk] = msloc;
    }
    //printf("\n");
    
    for(kk=0;kk<nholes;++kk){
      if(idh[kk] == -1){
        checkIonic = 1;
      }
    }

    //printf(" CheckIonic = %d\n",checkIonic);
    if(checkIonic == 0) {

      // Prepare address in base 6
      int addbase10,addbase6;
      addbase6 = 0;
      addbase10 = 0;
      for(kk=0;kk<nholes;++kk){
        //printf(" %d ",ms[kk]);
        addbase10 += ms[kk]*(int)pow(6,kk);
      }
      //printf("\n");
      addbase6 = base10ToBase6(addbase10);

      int idhole=0;
      for(kk=0;kk<nholes;++kk){
        idhole += idh[kk]*((int)pow(3,kk));
      }
      int idspin = charPtrToInt(dictionary_get(dadr,intToCharPtr(addbase6)));
      double fhole  = fachole[idhole];
      double fspin  = 1.0/charPtrToInt(dictionary_get(d, intToCharPtr(addbase6)));
      ix[0] = (size_t)iiii;
      //PetscCall(VecGetValues(valxr, 1, ix, wfy));
      VecGetValues(valxr, 1, ix, wfy);
      projvec[(iroot)*idx + idspin] += wfy[0]*fhole*sqrt(fspin);
      //printf(" %10.5f %10.5f\n",fhole*sqrt(fspin),wfy[0]);
      //*normproj += wfy[0]*wfy[0];
    }
    else{
      ix[0] = (size_t)iiii;
      //PetscCall(VecGetValues(valxr, 1, ix, y));
      VecGetValues(valxr, 1, ix, wfy);
      *normproj += wfy[0]*wfy[0];
      //printf(" y(%d)=%10.5f\n",iiii,wfy[0]);
    }
  }

  //printf(" norm=%10.5f\n",normproj[0]);
  //for(int ii=0;ii<idx;++ii){
  //  printf(" %10.5f\n",projvec[(iroot)*idx + ii]);
  //}

  if(nholes == 2){
    projvec2[(iroot)*idx + 0] = W_2holes;
  }
  else if(nholes == 3){
    projvec2[(iroot)*idx + 0] = W_3holes;
  }
  else if(nholes == 4){
    projvec2[(iroot)*idx + 0] = W_4holes;
  }
  else if(nholes == 5){
    projvec2[(iroot)*idx + 0] = W_5holes;
  }

  //printf("idxprojm=%d XS=%6.4f normproj=%6.4f \n",idxprojm,XS,normproj);
  //printf("--------- Hole Weight =%6.6f \n",projvec2[(iroot)*idx + 0]);
  //double sprojvec3=0;
  //double sprojvec1=0;
  //double sprojvec2=0;
  //for(int i=0;i<colm;++i){
  //  //sprojvec1 += projvec[(iroot)*idx + i];
  //  //sprojvec2 += projvec2[(iroot)*idx + i];
  //  sprojvec1 += projvec3[i];
  //  sprojvec2 += projvec4[i];
  //  //projvec2[(iroot)*idx + i] = fabs(projvec3[i])*fabs(projvec3[i]) + fabs(projvec4[i])*fabs(projvec4[i]);
  //  //printf(" (%8.4f %8.4f ) ",projvec[(iroot)*idx + i],projmatrix[idxprojm * colm + i]);
  //}
  //printf("\n");
  //printf("projvec1=%6.4f projvec2=%6.4f \n",sprojvec1,sprojvec2);
  //for(int i=0;i<colm;++i){
  //  printf(" (%8.4f %8.4f ) ",projvec2[(iroot)*idx + i],projmatrix[(idxprojm-1) * colm + i]);
  //}
  //printf("\n");
  //printf("proj total =%6.4f \n",sprojvec1*sprojvec1 + sprojvec2*sprojvec2);

  // Free dictionary
  dictionary_free(d);
  dictionary_free(dadr);

} /** END **/

