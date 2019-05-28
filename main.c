#include <stdio.h>
#include <stdlib.h>
//Ni uzas na mpq_t, de gmp, por kalkuli numbrojn realajn!
//Por ĉi biblioteko, uzu la argumenton -lgmp dum kiam konstruas
#include <gmp.h>

//Biblioteko por dum kiam ni skribas dosieron de png
//Por ĉi biblioteko, uzu la argumenton -lm dum kiam konstruas
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include "strukturoj.h"

#include "vektorfunkcioj_f.h"

void kalkulu_cxielkoloron(kolorof fk, duonrektoqq dr, vektoro_reala vktmp, mpf_t longo);
void kalkulu_cxielkoloron_f(kolorof fk, duonrektoqf *dr, vektorof vktmp);

void stbi_write_function(void *context, void *data, int size);

void kalkuluKoloron_duonrektoqq (kolorof fk, duonrektoqq dr, struct konstantoj *konstantoj);

void kreuBildonKolorprovanf (int *nx, int *ny, kolorof **bildof);

inline void kreuBildon255 (int *nx, int *ny, kolorof **bildof, koloro255 *bildo255);

void kreuBildonf (int *nx, int *ny, kolorof **bildof, char *reala);
void kreuBildonf_mpq (int *nx, int *ny, kolorof **bildof);
inline void kreuBildonf_mpz (int *nx, int *ny, kolorof **bildof, const unsigned long int pdm);

inline void skribuPPM (int *nx, int *ny, kolorof **bildof, FILE * eliro);

void skribuBildon (int *nx, int *ny, kolorof **bildof, FILE * eliro);

int main() {
  //larĝo
  static int nx = 400;
  //static int nx = 1280;
  //longo
  static int ny = 200;
  //static int ny = 720;
  char reala = 1; 

  //subtavolo de la bildo
  //kolorof bildof[ny][nx];
  kolorof **bildof = malloc(sizeof(*bildof)*ny);
  for (int y=0; y<ny; y++){
    bildof[y] = malloc(sizeof(*bildof[y])*nx);
  }
  //kreu la bildon
  kreuBildonf(&nx, &ny, bildof, &reala);
  //kreuBildonKolorprovanf(&nx, &ny, bildof);

  //skribu la bildon al stdout
  skribuBildon(&nx, &ny, bildof, stdout);

  //eliru
  for (int y=0; y<ny; y++){
    free(bildof[y]);
  }
  free(bildof);
  exit(EXIT_SUCCESS);
}


void kreuBildonf_mpq (int *nx, int *ny, kolorof **bildof){
  //kreu konstantojn
  struct konstantoj konstantoj;
  mpq_init(konstantoj.unu);
  mpq_init(konstantoj.du);
  mpq_init(konstantoj.kvar);
  mpq_set_ui(konstantoj.unu, 1, 1);
  mpq_set_ui(konstantoj.du, 2, 1);
  mpq_set_ui(konstantoj.kvar, 4, 1);
  
  int y, x, dux, duy;//, nny, nnx, dunnx, dunny, centnx, centny;
  //horizontala = (*nx)/100;
  //virtikala = (*ny)/100;
  //minimumau = -(horizontala/2);
  //minimumav = -(virtikala/2);
  
  kolorof nk;

  duonrektoqq ndr;
  duonrektonqq_init(ndr);

  //mpq_t longo;
  //mpq_init(longo);
  //dunx = (*nx)*2;
  //duny = (*ny)*2;
  //nnx = -(*nx);
  //nny = -(*ny);
  //dunnx = nnx*2;
  //dunny = nny*2;
  //centnx = (*nx)*100;
  //centny = (*ny)*100;

  for (y = 1; y <= *ny; y++){
    //v = y/(*nx);
    duy = y*2;
    for (x = 1; x <= *nx; x++){
      //u = x/(*nx);
      dux = x*2;
      //vektoron_realan_algxustigu(ndr[1], minimumau+(u*horizontala), minimumav+(v*virtikala), -1);
      //                            (-(horizontala/2))+((x/(*nx))*((*nx)/100))
      //                            ((x*(*nx))/(100*(*nx)))-(horizontala/2)
      //                            ((x*(*nx))/(100*(*nx)))-(((*nx)/100)/2)
      //                            ((x*(*nx))/(100*(*nx)))-((*nx)/200)
      //                            (x/100)-((*nx)/200)
      //                            (dux/200)-((*nx)/200)
      //                            (dux-(*nx))/200
      //
      //vektoron_realan_algxustigu_si(ndr[1], (nnx+dux), *nx, (nny+duy), *ny, -1, 1);
      vektoron_realan_algxustigu_si(ndr[1], (dux-(*nx)), 200, -(duy-(*ny)), 200, -1, 1);
      
      kalkuluKoloron_duonrektoqq(nk, ndr, &konstantoj);

      //if (y<=(*ny)/4){
      //	fprintf(stderr, "%f %f %f\n", nk[0], nk[1], nk[2]);
      //}
      
      bildof[y-1][x-1][0] = nk[0];
      bildof[y-1][x-1][1] = nk[1];
      bildof[y-1][x-1][2] = nk[2];
    }
  }
  duonrektonqq_klarigu(ndr);
  //mpq_clear(longo);
  mpq_clear(konstantoj.unu);
  mpq_clear(konstantoj.du);
  mpq_clear(konstantoj.kvar);
}

void kreuBildonf (int *nx, int *ny, kolorof **bildof, char *reala){
  //if (*reala == 0){
  //  kreuBildonf(nx, ny, bildof);
  //} else {
    kreuBildonf_mpq(nx, ny, bildof);
    //}
}

int kreuBildonf_mpz (int *nx, int *ny, kolorof **bildof, const unsigned long int pdm){
  //(pdm = potenco de metro, la potenco de 2 kiu egalas metron)
  int y, x;
  const mpz_t metro;
  mpz_ui_pow_ui(metro, 2, pdm);

  kolorof nk;
  duonrekto ndr;
  duonrekton_init(ndr);

  const mpz_t rastrumerlongo;
  
  
  for (y = 1; y <= *ny; y++){
    //v = y/(*nx);
    for (x = 1; x <= *nx; x++){
      //u = x/(*nx);
      
      //vektoron_algxustigu(ndr[1], minimumau+(u*horizontala), minimumav+(v*virtikala), -1);

    }
  }
}

/* void kalkulu_kolisiojn(kolisio *finakolisio, duonrektoqq dr, struct tmpvariabloj *tv, struct objektlisto *objektoj){ */
/*   //kalkulu kolisiojn kun objektaroj kaj kalkulu kolisiojn kun objektoj   */
/*   //kalkulu kolisiojn kun objekteroj */

/*   unsigned int i, j; */
/*   kolisioj kolisiojtmp; */
/*   kolisiotmp.numbro_da_kolisiovektoroj = 0; */
/*   //duonrektonqq_init(kolisiojtmp.kolisio.vektoro); */

/*   finakolisio -> estas = 0;  */
  
/*   //ĉu mi devus uzi x-on, y-on, aŭ z-on por kontroli la distancon? */
/*   int kontroldirekton, sgn, cmp; */
/*   mpq_abs(tv->vkrtmp[2], dr[1][2]); */
/*   mpq_abs(tv->vkrtmp[1], dr[1][1]); */
/*   mpq_abs(tv->vkrtmp[0], dr[1][0]); */
/*   if (mpq_cmp(tv->vkrtmp[2], tv->vkrtmp[1])>=0){ */
/*     //z estas pli granda ol y aŭ egalas ĝin */
/*     if (mpq_cmp(tv->vkrtmp[2], tv->vkrtmp[0])>=0){ */
/*       //z estas la plej granda! */
/*       kontroldirekton = 2; */
/*     } else { */
/*       //x estas pli granda ol z! Sed ĉu ĝi estas pli granda ol y? Hmm... */
/*       //Nu, y estas pli malgranda ol aŭ egalas z, tial x estas pli granda ol ĝi ankaŭ */
/*       kontroldirekton = 0; */
/*     } */
/*   } else { */
/*     //y estas pli granda ol z, sed ĉu ĝi estas pli granda ol x? Hmm... */
/*     if (mpq_cmp(tv->vkrtmp[1], tv->vkrtmp[0])>=0){ */
/*       //y estas pli granda aŭ egalas x kaj z */
/*       kontroldirekton = 1; */
/*     } else { */
/*       //x estas pli granda ol y kiu estas pli granda ol z */
/*       kontroldirekton = 0; */
/*     } */
/*   } */
/*   //nun ni havas ĝin! */
/*   //kalkulu se negativa */
/*   sgn = mpq_sgn(dr[1][kontroldirekton]); */
/*   //kalkulu kolisiojn en objektaroj */

/*   //unue, la objekteroj sferaj */
/*   for (i = 0; i < objektoj->numbro_da_objektaroj_sferaj; i++){ */
/*     //se ni havas kolision kun ĝi, kalkulu se ĝi estas pli proksima al nia kolisio */
/*     kalkulu_kolisiojn(&kolisiojtmp, dr, objektoj->objektaro_sfera[i].objektlisto); */
/*     if (kolisiotmp.estas > 0){ */
/*       if (sgn>=0){ */
/* 	//Konservu se estas pli granda, ĉar nia signo estas positiva */
/* 	cmp = mpq_cmp(kolisioj.kolisiejo[kontroldirekton], kolisiojtmp.kolisiejo[kontroldirekton]); */
/* 	if (cmp == 0){ */
/* 	  kolisioj->numbro_da_kolisiovektoroj += 1;  */
/* 	  if (kolisioj->numbro_da_kolisiovektoroj==0){ */
/* 	    kolisioj->kolisio = malloc(sizeof(kolisiovektoro)); */
/* 	    vektoron_realan_init(kolisioj->kolisio[0].vektoro); */
/* 	  } else { */
/* 	    realloc(kolisioj->kolisio, sizeof(kolisiovektoro)*(kolisioj->numbro_da_kolisiovektoroj)); */
/* 	    vektoron_realan_init(kolisioj->kolisio[(kolisioj->numbro_da_vektoroj)-1]); */
/* 	  } */
	  
/* 	} else if (cmp > 0) { */
/* 	  if (kolisioj->numbro_da_kolisiovektoroj>0){ */
/* 	    //forigu ilin – ni ne bezonas ilin */
/* 	    kolision_rekomencigu(kolisioj); */
/* 	  } */
/* 	  //ŝanĝŭ al la kolisio */
	  
/* 	}//else, ne estas grava */
/*       } else if (sgn<=0) { */
/* 	//Konservu se estas malpli granda, ĉar nia signo estas negativa */
/*       } */
/*     } */
/*   } */

/*   duonrektonqq_klarigu(kolisiojtmp.kolisio.vektoro); */
/* } */

void kalkuluKoloron_duonrektoqq(kolorof fk, duonrektoqq dr, struct konstantoj *konstantoj){
  //fk estas la Fina Koloro


  //kreu variablojn dumtempajn
  struct tmpvariabloj tv;
  tmpvariablojn_init(&tv);
  //struct tmpvariabloj {
  //  mpq_t mpqtmp;
  //  mpq_t mpqtmp2;
  //  mpf_t mpftmp;
  //  vektoro_reala vkrtmp;
  //  vektoro_reala vkrtmp2;
  //};
  
  vektoro_reala lmc, centrodesfero;
  vektoron_realan_init(lmc);
  vektoron_realan_init(centrodesfero);
  duonrektoqf dr_qf;
  dr_qf.loko = &(dr[0]);
  dr_qf.vkt[0] = (float) mpq_get_d(dr[1][0]);
  dr_qf.vkt[1] = (float) mpq_get_d(dr[1][1]);
  dr_qf.vkt[2] = (float) mpq_get_d(dr[1][2]);
  vektorof vktf;
  vektorompf vktmpf;
  vektoronmpf_init(vktmpf);
  mpq_t radiuso;
  //mpqtmp2 = &(tv.mpqtmp2);
  mpq_init(radiuso);
  mpq_set_ui(radiuso, 1, 2);
  mpf_t respondo;
  mpf_init(respondo);


  mpq_set_si(centrodesfero[2], -1, 1);

  //kalkulu kolisiojn
  //kalkulu_kolisiojn(mondo);
  duonrektoqq_tusxas_sferon(respondo, centrodesfero, radiuso, dr, lmc, &tv, konstantoj);

  if (mpf_cmp_d(respondo, 0.0)>0){
    //montru_al_loko(duonrektoqq dr, t) { return dr[0] + t*dr[1]; }
    //vktmp = vektoro_de_unuo(montru_al_loko(dr, respondo)-vkt{0, 0 -1));
    mpq_set_f(tv.mpqtmp, respondo);
    montru_al_loko_qq(tv.vkrtmp, dr, tv.mpqtmp);
    //mpq_add_si(vktmp[2], vktmp[2], konstantoj->unu);
    //vktmp = vktmp - vkt{0, 0, -1}
    mpz_add(mpq_numref(tv.vkrtmp[2]), mpq_numref(tv.vkrtmp[2]), mpq_denref(tv.vkrtmp[2]));
    //vektoro_de_unuo
    longo_de_vektoro_reala_vktmp(tv.mpftmp, tv.vkrtmp, tv.vkrtmp2);
    
    vektorompq_al_vektorompf(vktmpf, tv.vkrtmp);
    vektoronmpf_dividu_per_mpf(vktmpf, vktmpf, tv.mpftmp);
    //vktmpf nun estas eta
    vektorompf_al_vektorof(vktf, vktmpf);
    
    //kalkulu per eta vktmpf
    vektorof_mul(vktf, vktf, 0.5);
    vektorof_aldonu(vktf, vktf, 0.5);
    
    //fino
    fk[0] = vktf[0];
    fk[1] = vktf[1];
    fk[2] = vktf[2];

    //(malnova)
    //mpq_set_f(mpqtmp2, mpftmp);
    //mpq_inv(mpqtmp2, mpqtmp2);
    //vektoron_realan_mul_per_numbro(vktmp, vktmp, mpqtmp2);
    //
    //
    ////return 0.5*(vktmp+1)
    //mpz_add(mpq_numref(vktmp[0]), mpq_numref(vktmp[0]), mpq_denref(vktmp[0]));
    //mpz_add(mpq_numref(vktmp[1]), mpq_numref(vktmp[1]), mpq_denref(vktmp[1]));
    //mpz_add(mpq_numref(vktmp[2]), mpq_numref(vktmp[2]), mpq_denref(vktmp[2]));
    //
    //mpz_add(mpq_denref(vktmp[0]), mpq_denref(vktmp[0]), mpq_denref(vktmp[0]));
    //mpz_add(mpq_denref(vktmp[1]), mpq_denref(vktmp[1]), mpq_denref(vktmp[1]));
    //mpz_add(mpq_denref(vktmp[2]), mpq_denref(vktmp[2]), mpq_denref(vktmp[2]));
    //
    //
    //fk[0] = (float) mpq_get_d(vktmp[0]);
    //fk[1] = (float) mpq_get_d(vktmp[1]);
    //fk[2] = (float) mpq_get_d(vktmp[2]);
  } else {
    //kalkulu_cxielkoloron(fk, dr, vktmp, respondo);
    kalkulu_cxielkoloron_f(fk, &dr_qf, vktf);
  }

  
  vektoron_realan_klarigu(lmc);
  vektoron_realan_klarigu(centrodesfero);
  vektoronmpf_klarigu(vktmpf);
  mpq_clear(radiuso);
  mpf_clear(respondo);
  tmpvariablojn_klarigu(&tv);
}

void kalkulu_cxielkoloron_f(kolorof fk, duonrektoqf *dr, vektorof vktmp){
      
  //mi uzas matimatikon

  //y = dr[1][1]
  //longo_de_vektoro(dr) = sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2)
  //vektoro_de_unuo(dr).y = y / longo_de_vektoro(dr)
  //vektoro_de_unuo(dr).y = y / sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2)
  //ud.y = vektoro_de_unuo(dr).y
  //ud.y = y / sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2)
  //t = 0.5*(ud.y + 1)
  //t = 0.5*((y / sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2) + 1)
  //
  //(1.0-t)*vektoro(0.5, 0.7, 1.0) + t*vektoro(1, 1, 1)
  //
  //

  //nun, ni kalkulu!
  //variablojn dumtempajn
    //mpf_t longo, dulongo, longopy, longomy, tmpf;
  float longo, dulongo, longopy, longomy;//, tmpf;
  //mpf_set_default_prec(32);
  //mpf_init2(longo, 64);
  //mpf_init2(dulongo, 64);
  //mpf_init2(longopy, 64);
  //mpf_init2(longomy, 64);
  //mpf_init2(tmpf, 64);
  
  //sqrt de x^2+y^2+z^2 -- ni uzu na longo_de_vektoro_vktmp(mpf_t longo, vektoro_reala dr[1] "vektoro", vektoro_reala vktmp "dumtempa vektoro")
  longo_de_vektorof_vktmp(&longo, dr->vkt, vktmp);
  //fprintf(stderr, "%f\n", dr->vkt[1]);

  //konvertu na y al mpf en tmpf
  //mpf_set_q(tmpf, dr[1][1]);
  
  
  //t = 0.5*((y/longo)+1)
  //3 op
  //t = 0.5*((y+longo)/longo)
  //t = (y+longo)/(2*longo)
  //3 op
  //(1.0-t)*vektoro(0.5, 0.7, 1.0) + t*vektoro(1, 1, 1)
  //7 op
  //(1.0-((y+longo)/(2*longo)))*vektoro(0.5, 0.7, 1.0) +
  //((y+longo)/(2*longo))*vektoro(0.5, 0.7, 1.0)
  //
  //( ((-y-longo)/(2*longo)) + (1/1) ) * vektoro(0.5, 0.7, 1.0)
  //( ((-y-longo)/(2*longo)) + (2*longo/2*longo) ) * vektoro(0.5, 0.7, 1.0)
  //( ((2*longo)-y-longo) / (2*longo) ) * vektoro(0.5, 0.7, 1.0)
  //( ((longo+longo)-y-longo) / (2*longo) ) * vektoro(0.5, 0.7, 1.0)
  //( (longo-y)/(2*longo) ) * vektoro(0.5, 0.7, 1.0) +
  //( (longo+y)/(2*longo) ) * vektoro(1, 1, 1)
  //8 op each
  //
  //fk[0] = ( (0.5*(longo-y))/(2*longo) ) +
  //        ( (longo+y)/(2*longo) )
  // OR
  //fk[0] = ( (longo-y)/(4*longo) ) +
  //        ( (longo+y)/(2*longo) )
  //mi ŝatas la unuan
  //fk[0] = ( (longo-y)+0.5*(longo+y) )/(2*longo)

  //kalkulu na longo-y
  //mpf_sub(longomy, longo, tmpf);
  longomy = longo - (dr->vkt[1]);
  //kalkulu na longo+y
  //mpf_add(longopy, longo, dr[1][1]);
  longopy = longo + (dr->vkt[1]);
  //kalkulu na dulongo
  //mpf_mul_ui(dulongo, longo, 2);
  dulongo = longo*2;

  //nun, kalkulu na fk[0]
  //unue, alĝustigu tmpf al 0.5
  //mpf_set_d(tmpf, 0.5);
  //mpf_mul(tmpf, longomy, tmpf);
  fk[0] = longopy*0.5;
  //mpf_add(tmpf, tmpf, longopy);
  fk[0] = fk[0] + longomy;
  //mpf_div(tmpf, tmpf, dulongo);
  fk[0] = fk[0]/dulongo;
  //fk[0] = mpf_get_d(tmpf);
  
  //fk[1] = ( (0.7*(longo-y))/(2*longo) ) +
  //        ( (longo+y)/(2*longo) )
  //fk[1] = ( (0.7*(longo-y)+(longo+y))/(2*longo)

  //kalkulu fk[1]
  //unue, alĝustigu tmpf al 0.7
  //mpf_set_d(tmpf, 0.7);
  //mpf_mul(tmpf, longomy, tmpf);
  fk[1] = longopy*0.7;
  //gmp_fprintf(stderr, "%s\n", mpf_get_str(NULL, &exp2, 10, 5, tmpf));
  //mpf_add(tmpf, tmpf, longopy);
  fk[1] = fk[1] + longomy;
  //mpf_div(tmpf, tmpf, dulongo);
  fk[1] = fk[1] / dulongo;
  //fk[1] = mpf_get_d(tmpf);
  
  //fk[2] = ( (longo-y)/(2*longo) ) 
  //        ( (longo+y)/(2*longo) )
  //fk[2] = ( (2*longo)/(2*longo) )
  //fk[2] = 1
  fk[2] = 1;

  //klarigu
  //mpf_clear(longo);
  //mpf_clear(dulongo);
  //mpf_clear(longopy);
  //mpf_clear(longomy);
  
  //mpf_clear(tmpf);
}

void kalkulu_cxielkoloron(kolorof fk, duonrektoqq dr, vektoro_reala vktmp, mpf_t longo){
    
  //mi uzas matimatikon

  //y = dr[1][1]
  //longo_de_vektoro(dr) = sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2)
  //vektoro_de_unuo(dr).y = y / longo_de_vektoro(dr)
  //vektoro_de_unuo(dr).y = y / sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2)
  //ud.y = vektoro_de_unuo(dr).y
  //ud.y = y / sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2)
  //t = 0.5*(ud.y + 1)
  //t = 0.5*((y / sqrt(dr[0]^2 + dr[1]^2 + dr[2]^2) + 1)
  //
  //(1.0-t)*vektoro(1, 1, 1) + t*vektoro(0.5, 0.7, 1.0)
  //
  //

  //nun, ni kalkulu!
  //variablojn dumtempajn
  mpf_t dulongo, longopy, longomy, tmpf;
  mpf_set_default_prec(32);
  //mpf_init2(longo, 64);
  mpf_init2(dulongo, 64);
  mpf_init2(longopy, 64);
  mpf_init2(longomy, 64);
  mpf_init2(tmpf, 64);
  
  //sqrt de x^2+y^2+z^2 -- ni uzu na longo_de_vektoro_vktmp(mpf_t longo, vektoro_reala dr[1] "vektoro", vektoro_reala vktmp "dumtempa vektoro")
  longo_de_vektoro_reala_vktmp(longo, dr[1], vktmp);

  //konvertu na y al mpf en tmpf
  mpf_set_q(tmpf, dr[1][1]);
  
  
  //t = 0.5*((y/longo)+1)
  //3 op
  //t = 0.5*((y+longo)/longo)
  //t = (y+longo)/(2*longo)
  //3 op
  //(1.0-t)*vektoro(0.5, 0.7, 1.0) + t*vektoro(1, 1, 1)
  //7 op
  //(1.0-((y+longo)/(2*longo)))*vektoro(0.5, 0.7, 1.0) +
  //((y+longo)/(2*longo))*vektoro(0.5, 0.7, 1.0)
  //
  //( ((-y-longo)/(2*longo)) + (1/1) ) * vektoro(0.5, 0.7, 1.0)
  //( ((-y-longo)/(2*longo)) + (2*longo/2*longo) ) * vektoro(0.5, 0.7, 1.0)
  //( ((2*longo)-y-longo) / (2*longo) ) * vektoro(0.5, 0.7, 1.0)
  //( ((longo+longo)-y-longo) / (2*longo) ) * vektoro(0.5, 0.7, 1.0)
  //( (longo-y)/(2*longo) ) * vektoro(0.5, 0.7, 1.0) +
  //( (longo+y)/(2*longo) ) * vektoro(1, 1, 1)
  //8 op each
  //
  //fk[0] = ( (0.5*(longo-y))/(2*longo) ) +
  //        ( (longo+y)/(2*longo) )
  // OR
  //fk[0] = ( (longo-y)/(4*longo) ) +
  //        ( (longo+y)/(2*longo) )
  //mi ŝatas la unuan
  //fk[0] = ( 0.5*(longo-y)+(longo+y) )/(2*longo)

  //kalkulu na longo-y
  mpf_sub(longomy, longo, tmpf);
  //kalkulu na longo+y
  mpf_add(longopy, longo, tmpf);
  //kalkulu na dulongo
  mpf_mul_ui(dulongo, longo, 2);

  //nun, kalkulu na fk[0]
  //unue, alĝustigu tmpf al 0.5
  mpf_set_d(tmpf, 0.5);
  mpf_mul(tmpf, longomy, tmpf);
  mpf_add(tmpf, tmpf, longopy);
  mpf_div(tmpf, tmpf, dulongo);
  fk[0] = mpf_get_d(tmpf);
  
  //fk[1] = ( (0.7*(longo-y))/(2*longo) ) +
  //        ( (longo+y)/(2*longo) )
  //fk[1] = ( (0.7*(longo-y)+(longo+y))/(2*longo)

  //kalkulu fk[1]
  //unue, alĝustigu tmpf al 0.7
  mpf_set_d(tmpf, 0.7);
  mpf_mul(tmpf, longomy, tmpf);
  //gmp_fprintf(stderr, "%s\n", mpf_get_str(NULL, &exp2, 10, 5, tmpf));
  mpf_add(tmpf, tmpf, longopy);
  mpf_div(tmpf, tmpf, dulongo);
  fk[1] = mpf_get_d(tmpf);
  
  //fk[2] = ( (longo-y)/(2*longo) ) 
  //        ( (longo+y)/(2*longo) )
  //fk[2] = ( (2*longo)/(2*longo) )
  //fk[2] = 1
  fk[2] = 1;

  //klarigu
  //mpf_clear(longo);
  mpf_clear(dulongo);
  mpf_clear(longopy);
  mpf_clear(longomy);
  mpf_clear(tmpf);
}

void kreuBildonKolorprovanf (int *nx, int *ny, kolorof **bildof){
  int y;
  int x;
  for (y = 0; y < *ny; y++){
    for (x = 0; x < *nx; x++){
      bildof[y][x][0] = ((float)x) / ((float)((*nx)-1));
      bildof[y][x][1] = ((float)(((*ny)-1)-y))/((float)((*ny)-1));
      bildof[y][x][2] = 0.2;
    }
  }
}

void skribuBildon (int *nx, int *ny, kolorof **bildof, FILE * eliro){
  //skribuPPM(nx, ny, bildof, eliro);
  //koloro255 bildo255[*ny][*nx];
  koloro255 *bildo255 = malloc(sizeof(*bildo255)*(*ny)*(*nx));
  kreuBildon255(nx, ny, bildof, bildo255);
  stbi_write_png_to_func(&stbi_write_function, eliro, *nx, *ny, 3, bildo255, sizeof(bildo255[0][0])*(*nx)*3);
  //for (int y=0; y<(*ny); y++){
    //free(bildo255[y]);
    //}
  free(bildo255);
}

void kreuBildon255 (int *nx, int *ny, kolorof **bildof, koloro255 *bildo255){
  int y, x;
  for (y = 0; y < (*ny); y++){
    for (x = 0; x < (*nx); x++){
      bildo255[(y*(*nx))+x][0] = (unsigned char)(255*bildof[y][x][0]);
      bildo255[(y*(*nx))+x][1] = (unsigned char)(255*bildof[y][x][1]);
      bildo255[(y*(*nx))+x][2] = (unsigned char)(255*bildof[y][x][2]);
    }
  }
}

void skribuPPM (int *nx, int *ny, kolorof **bildof, FILE * eliro){
  fprintf(eliro, "P3\n%d %d 255\n", *nx, *ny);
  koloro255 kol;
  for (int j = 0; j < *ny; j++){
    for (int i = 0; i < *nx; i++){
      kol[0] = (unsigned char)(255*bildof[j][i][0]);
      kol[1] = (unsigned char)(255*bildof[j][i][1]);
      kol[2] = (unsigned char)(255*bildof[j][i][2]);
      fprintf(eliro, "%hhu %hhu, %hhu\n", kol[0], kol[1], kol[2]);
    }
  }
}

void stbi_write_function(void *context, void *data, int size){
  fwrite(data, 1, size, context);
}
