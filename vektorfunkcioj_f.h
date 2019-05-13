#include "strukturoj.h"

#ifndef _VEKTORFUNKCIOJ_F_H_
#define _VEKTORFUNKCIOJ_F_H_


//inline void vektoronf_init(vektorof vkt);
inline void vektoronmpf_init(vektorompf vkt);
inline void duonrektonqf_init(duonrektoqf *dr);

//inline void vektoronf_klarigu(vektorof vkt);
inline void vektoronmpf_klarigu(vektorompf vkt);
inline void duonrektonqf_klarigu(duonrektoqf *dr);

inline void vektoronf_algxustigu(vektorof vkt, float x, float y, float z);

inline void vektorompq_al_vektorompf(vektorompf vktmpf, vektoro_reala vkt);

inline void vektorompf_al_vektorof(vektorof vktf, vektorompf vktmpf);

inline void vektoronmpf_dividu_per_mpf(vektorompf vkt, vektorompf vkten, mpf_t dividanto);

inline void vektorof_aldonu(vektorof vktf, vektorof vktfen, float aldonanto);

inline void vektorof_mul(vektorof vktf, vektorof vktfen, float multiplikanto);

inline void longo_de_vektorof_vktmp(float *longo, vektorof vkt, vektorof vktmp);

inline void vektoronf_mul_per_vektoronumbrojf(vektorof vkt_fina, const vektorof vkt1, const vektorof vkt2);

void vektoronf_dot(float *dot, vektorof vkt1, vektorof vkt2);

//int duonrektoqf_tusxas_sferon(vektoro centrodesfero, mpq_t radiuso, const duonrektoqf dr);

void kalkuluKoloron_duonrektoqf (kolorof fk, duonrektoqf *dr);

inline void kreuBildonf_f (int *nx, int *ny, kolorof **bildof);

void kalkuluKoloron_duonrektoqf (kolorof fk, duonrektoqf *dr){
  //fk estas la Fina Koloro

  
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
  //kreu variablojn dumtempajn

  //por mpf_t
  vektorof vktmp; 
  //vektoronf_init(vktmp);
  //mpf_t longo, dulongo, longopy, longomy, tmpf;
  float longo, dulongo, longopy, longomy;//, tmpf;
  //mpf_set_default_prec(32);
  //mpf_init2(longo, 64);
  //mpf_init2(dulongo, 64);
  //mpf_init2(longopy, 64);
  //mpf_init2(longomy, 64);
  //mpf_init2(tmpf, 64);
  
  //sqrt de x^2+y^2+z^2 -- ni uzu na longo_de_vektoro_vktmp(mpf_t longo, vektoro dr[1] "vektoro", vektoro vktmp "dumtempa vektoro")
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
  //fk[0] = ( 0.5*(longo-y)+(longo+y) )/(2*longo)

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
  fk[0] = longomy*0.5;
  //mpf_add(tmpf, tmpf, longopy);
  fk[0] = fk[0] + longopy;
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
  fk[1] = longomy*0.7;
  //gmp_fprintf(stderr, "%s\n", mpf_get_str(NULL, &exp2, 10, 5, tmpf));
  //mpf_add(tmpf, tmpf, longopy);
  fk[1] = fk[1] + longopy;
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
  //vektoronf_klarigu(vktmp);

}

void vektoronmpf_dividu_per_mpf(vektorompf vkt, vektorompf vkten, mpf_t dividanto){
  mpf_div(vkt[0], vkten[0], dividanto);
  mpf_div(vkt[1], vkten[1], dividanto);
  mpf_div(vkt[2], vkten[2], dividanto);
}

void vektorof_mul(vektorof vktf, vektorof vktfen, float multiplikanto){
  vktf[0] = vktfen[0] * multiplikanto;
  vktf[1] = vktfen[1] * multiplikanto;
  vktf[2] = vktfen[2] * multiplikanto;
}

void vektorof_aldonu(vektorof vktf, vektorof vktfen, float aldonanto){
  vktf[0] = vktfen[0] + aldonanto;
  vktf[1] = vktfen[1] + aldonanto;
  vktf[2] = vktfen[2] + aldonanto;
}

void vektorompf_al_vektorof(vektorof vktf, vektorompf vktmpf){
  vktf[0] = mpf_get_d(vktmpf[0]);
  vktf[1] = mpf_get_d(vktmpf[1]);
  vktf[2] = mpf_get_d(vktmpf[2]);
}

void vektorompq_al_vektorompf(vektorompf vktmpf, vektoro_reala vkt){
  mpf_set_q(vktmpf[0], vkt[0]);
  mpf_set_q(vktmpf[1], vkt[1]);
  mpf_set_q(vktmpf[2], vkt[2]);
}

void longo_de_vektorof_vktmp(float *longo, vektorof vkt, vektorof vktmp){
  vektoronf_mul_per_vektoronumbrojf(vktmp, vkt, vkt);
  //aldonu ĉiujn numbrojn en vktmp al vktmp[0]
  vktmp[0] = vktmp[0] + vktmp[1];
  vktmp[0] = vktmp[0] + vktmp[2];

  //sqrt longon
  *longo = sqrt(vktmp[0]);
}

void vektoronf_dot(float *dot, vektorof vkt1, vektorof vkt2){
  vektorof vktmp;
  vektoronf_mul_per_vektoronumbrojf(vktmp, vkt1, vkt2);
  (*dot) = vktmp[0] + vktmp[1];
  (*dot) = (*dot) + vktmp[2];
}

void vektoronf_mul_per_vektoronumbrojf(vektorof vkt_fina, const vektorof vkt1, const vektorof vkt2){
  vkt_fina[0] = vkt1[0]*vkt2[0];
  vkt_fina[1] = vkt1[1]*vkt2[1];
  vkt_fina[2] = vkt1[2]*vkt2[2];
}

void duonrektonqf_klarigu(duonrektoqf *dr){
  vektoron_realan_klarigu(*(dr->loko));
  //vektoronf_klarigu(dr->vkt);
}

void duonrektonqf_init(duonrektoqf *dr){
  vektoro_reala drloko;
  dr->loko = &drloko;
  vektoron_realan_init(*(dr->loko));
  //vektoronf_init(dr->vkt);
}

void vektoronmpf_init(vektorompf vkt){
  mpf_init(vkt[0]);
  mpf_init(vkt[1]);
  mpf_init(vkt[2]);
}
//void vektoronf_init(vektorof vkt){
  //faru nenion
//}

void vektoronmpf_klarigu(vektorompf vkt){
  mpf_clear(vkt[0]);
  mpf_clear(vkt[1]);
  mpf_clear(vkt[2]);
}
//void vektoronf_klarigu(vektorof vkt){
  //faru nenion
//}


void vektoronf_algxustigu(vektorof vkt, float x, float y, float z){
  vkt[0] = x;
  vkt[1] = y;
  vkt[2] = z;
}

void kreuBildonf_f (int *nx, int *ny, kolorof **bildof){
  int y, x;//, nny, nnx, dux, duy;
  float u, v, horizontala, virtikala, minimumau, minimumav;
  horizontala = (*nx)/100;
  virtikala = (*ny)/100;
  minimumau = -(horizontala/2);
  minimumav = -(virtikala/2);

  kolorof nk;
  
  duonrektoqf ndr;
  duonrektonqf_init(&ndr);

  
  //mpq_t longo;
  //mpq_init(longo);
  //dunx = (*nx)*2;
  //duny = (*ny)*2;
  //nnx = -(*nx);
  //nny = -(*ny);

  for (y = 1; y <= *ny; y++){
    //duy = y*2;
    v = (float)y/(float)(*ny);
    for (x = 1; x <= *nx; x++){
      //dux = x*2;
      u = (float)x/(float)(*nx);
      vektoronf_algxustigu(ndr.vkt, minimumau+(u*horizontala), minimumav+(v*virtikala), -1);
      //fprintf(stderr, "%f\n", v);
      
      kalkuluKoloron_duonrektoqf(nk, &ndr);

      //if (y<=(*ny)/4){
      //	fprintf(stderr, "%f %f %f\n", nk[0], nk[1], nk[2]);
      //}
      
      bildof[y-1][x-1][0] = nk[0];
      bildof[y-1][x-1][1] = nk[1];
      bildof[y-1][x-1][2] = nk[2];
    }
  }
  duonrektonqf_klarigu(&ndr);
  //mpq_clear(longo);
}

#endif
