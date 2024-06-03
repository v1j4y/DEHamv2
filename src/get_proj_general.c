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
  //double fachole1[9];
  //double fachole2[9];
  //double fachole3[9];
  //double fachole4[9];
  //double fachole5[9];
  //double fachole6[9];
  //double fachole7[9];
  //double fachole8[9];
  //double fachole9[9];
  //double normproj1tot=0.0;
  //double normproj2tot=0.0;
  //double normproj3tot=0.0;
  //double normproj4tot=0.0;
  //double normproj5tot=0.0;
  //double normproj6tot=0.0;
  //double normproj7tot=0.0;
  //double normproj8tot=0.0;
  //double normproj9tot=0.0;
  //double normproj1[6];
  //double normproj2[6];
  //double normproj3[6];
  //double normproj4[6];
  //double normproj5[6];
  //double normproj6[6];
  //double normproj7[6];
  //double normproj8[6];
  //double normproj9[6];
  //for(int i=0;i<6;++i) {
  //  normproj1[i]=0.0;
  //  normproj2[i]=0.0;
  //  normproj3[i]=0.0;
  //  normproj4[i]=0.0;
  //  normproj5[i]=0.0;
  //  normproj6[i]=0.0;
  //  normproj7[i]=0.0;
  //  normproj8[i]=0.0;
  //  normproj9[i]=0.0;
  //}

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
  // First vector
  //fachole1[0] = 1./4;
  //fachole1[1] = -sqrt(1./8);
  //fachole1[2] = 1./4;
  //fachole1[3] = -sqrt(1./8);
  //fachole1[4] = 1./2;
  //fachole1[5] = -sqrt(1./8);
  //fachole1[6] = 1./4;
  //fachole1[7] = -sqrt(1./8);
  //fachole1[8] = 1./4;
  //// 2 vector
  //fachole2[0] = 1./4;
  //fachole2[1] = sqrt(1./8);
  //fachole2[2] = 1./4;
  //fachole2[3] = -sqrt(1./8);
  //fachole2[4] = -1./2;
  //fachole2[5] = -sqrt(1./8);
  //fachole2[6] = 1./4;
  //fachole2[7] = sqrt(1./8);
  //fachole2[8] = 1./4;
  //// 3 vector
  //fachole3[0] = -sqrt(1./8);
  //fachole3[1] = 0.0;
  //fachole3[2] = sqrt(1./8);
  //fachole3[3] = 1./2;
  //fachole3[4] = 0.0;
  //fachole3[5] = -1./2;
  //fachole3[6] = -sqrt(1./8);
  //fachole3[7] = 0.0;
  //fachole3[8] = sqrt(1./8);
  //// 4 vector
  //fachole4[0] = 1./4;
  //fachole4[1] = -sqrt(1./8);
  //fachole4[2] = 1./4;
  //fachole4[3] = sqrt(1./8);
  //fachole4[4] = -1./2;
  //fachole4[5] = sqrt(1./8);
  //fachole4[6] = 1./4;
  //fachole4[7] = -sqrt(1./8);
  //fachole4[8] = 1./4;
  //// 5 vector
  //fachole5[0] = 1./4;
  //fachole5[1] = sqrt(1./8);
  //fachole5[2] = 1./4;
  //fachole5[3] = sqrt(1./8);
  //fachole5[4] = 1./2;
  //fachole5[5] = sqrt(1./8);
  //fachole5[6] = 1./4;
  //fachole5[7] = sqrt(1./8);
  //fachole5[8] = 1./4;
  //// 6 vector
  //fachole6[0] = -sqrt(1./8);
  //fachole6[1] = 0.0;
  //fachole6[2] = sqrt(1./8);
  //fachole6[3] = -1./2;
  //fachole6[4] = 0.0;
  //fachole6[5] = 1./2;
  //fachole6[6] = -sqrt(1./8);
  //fachole6[7] = 0.0;
  //fachole6[8] = sqrt(1./8);
  //// 7 vector
  //fachole7[0] = -sqrt(1./8);
  //fachole7[1] = 1./2;
  //fachole7[2] = -sqrt(1./8);
  //fachole7[3] = 0.0;
  //fachole7[4] = 0.0;
  //fachole7[5] = 0.0;
  //fachole7[6] = sqrt(1./8);
  //fachole7[7] = -1./2;
  //fachole7[8] = sqrt(1./8);
  //// 8 vector
  //fachole8[0] = -sqrt(1./8);
  //fachole8[1] = -1./2;
  //fachole8[2] = -sqrt(1./8);
  //fachole8[3] = 0.0;
  //fachole8[4] = 0.0;
  //fachole8[5] = 0.0;
  //fachole8[6] = sqrt(1./8);
  //fachole8[7] = 1./2;
  //fachole8[8] = sqrt(1./8);
  //// 9 vector
  //fachole9[0] = 1./2;
  //fachole9[1] = 0.0;
  //fachole9[2] = -1./2;
  //fachole9[3] = 0.0;
  //fachole9[4] = 0.0;
  //fachole9[5] = 0.0;
  //fachole9[6] = -1./2;
  //fachole9[7] = 0.0;
  //fachole9[8] = 1./2;

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
    for(kk=0;kk<nholes;++kk) {
      idh[kk]=-1;
    }
   
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
      //normproj1[idspin] += (wfy[0]*fachole1[idhole]*sqrt(fspin));
      //normproj2[idspin] += (wfy[0]*fachole2[idhole]*sqrt(fspin));
      //normproj3[idspin] += (wfy[0]*fachole3[idhole]*sqrt(fspin));
      //normproj4[idspin] += (wfy[0]*fachole4[idhole]*sqrt(fspin));
      //normproj5[idspin] += (wfy[0]*fachole5[idhole]*sqrt(fspin));
      //normproj6[idspin] += (wfy[0]*fachole6[idhole]*sqrt(fspin));
      //normproj7[idspin] += (wfy[0]*fachole7[idhole]*sqrt(fspin));
      //normproj8[idspin] += (wfy[0]*fachole8[idhole]*sqrt(fspin));
      //normproj9[idspin] += (wfy[0]*fachole9[idhole]*sqrt(fspin));
      //printf(" %10.5f %10.5f\n",fhole*sqrt(fspin),wfy[0]);
      //printf("%d 1 idhole=%d fac=%8f fac1=%8f\n",iii,idhole,fachole[idhole],fachole2[idhole]);
    }
    else{
      ix[0] = (size_t)iiii;
      //PetscCall(VecGetValues(valxr, 1, ix, y));
      VecGetValues(valxr, 1, ix, wfy);
      normproj[0] += wfy[0]*wfy[0];
      //printf(" y(%d)=%10.5f\n",iiii,wfy[0]);
      //printf("%d 0 \n",iii);
    }
  }

  //for(int i=0;i<6;++i) {
  //  normproj1tot += normproj1[i]*normproj1[i];
  //  normproj2tot += normproj2[i]*normproj2[i];
  //  normproj3tot += normproj3[i]*normproj3[i];
  //  normproj4tot += normproj4[i]*normproj4[i];
  //  normproj5tot += normproj5[i]*normproj5[i];
  //  normproj6tot += normproj6[i]*normproj6[i];
  //  normproj7tot += normproj7[i]*normproj7[i];
  //  normproj8tot += normproj8[i]*normproj8[i];
  //  normproj9tot += normproj9[i]*normproj9[i];
  //}

  //printf(" norm=%8f norm2=%8f norm3=%8f sum=%8f\n",normproj1tot,normproj2tot,normproj3tot,normproj1tot+normproj2tot+normproj3tot);
  //printf(" norm=%8f norm2=%8f norm3=%8f sum=%8f\n",normproj4tot,normproj5tot,normproj6tot,normproj4tot+normproj5tot+normproj6tot);
  //printf(" norm=%8f norm2=%8f norm3=%8f sum=%8f\n",normproj7tot,normproj8tot,normproj9tot,normproj7tot+normproj8tot+normproj9tot);
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

