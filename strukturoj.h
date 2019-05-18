#ifndef _STRUCTUROJ_H_
#define _STRUCTUROJ_H_
#include <gmp.h>
typedef mpz_t vektoro[3];
typedef mpq_t vektoro_reala[3];
typedef vektoro_reala duonrektoqq[2];
typedef struct materialo {
  char tipo[64];
  void *agordoj;
} materialo;
typedef struct kolisiovektoro {
  vektoro vektoro;
  materialo *materialo;
} kolisiovektoro;
typedef struct kolisioj {
  unsigned char estas;
  vektoro kolisiejo;
  size_t numbro_da_kolisiovektoroj;
  kolisiovektoro *kolisio;
} kolisio;
typedef unsigned char koloro255[3];
typedef float kolorof[3];
//typedef struct bildof {
//  int nx;
//  int ny;
//  koloro rastrumero[ny][nx];
//};
typedef float vektorof[3];

typedef mpf_t vektorompf[3];

typedef struct duonrektoqf {
  vektoro_reala *loko;
  vektorof vkt;
} duonrektoqf;

typedef struct duonrekto {
  vektoro loko;
  vektorof direkcio;
} duonrekto;

struct sfero_reala {
  vektoro_reala centrodesfero;
  mpq_t radiuso;
};

struct sfero {
  vektoro centro;
  mpz_t radiuso;
  materialo *materialo;
};

struct objektlisto {
  char sxargxa;
  
  unsigned int numbro_da_objektaroj_sferaj;
  struct objektaro_sfera *objektaro_sfera;
  
  unsigned int numbro_da_sferoj;
  struct sfero *sferoj;
};

struct objektaro_sfera {
  struct objektlisto *objektlisto;
  struct sfero *objektarejo;
};



struct konstantoj {
  mpq_t unu;
  mpq_t du;
  mpq_t kvar;
};

struct tmpvariabloj {
  mpq_t mpqtmp;
  mpq_t mpqtmp2;
  mpz_t mpztmp;
  mpz_t mpztmp2;
  mpf_t mpftmp;
  vektoro_reala vkrtmp;
  vektoro_reala vkrtmp2;
  vektoro vktmp;
  vektoro vktmp2;
};



inline void vektoron_init(vektoro vkt);
inline void vektoron_realan_init(vektoro_reala vkt);
inline void duonrektonqq_init(duonrektoqq dr);

inline void vektoron_klarigu(vektoro vkt);
inline void vektoron_realan_klarigu(vektoro_reala vkt);
inline void duonrektonqq_klarigu(duonrektoqq dr);

inline void vektoron_algxustigu(vektoro vkt, mpz_t x, mpz_t y, mpz_t z);
inline void vektoron_realan_algxustigu(vektoro_reala vkt, mpq_t x, mpq_t y, mpq_t z);
inline void vektoron_algxustigu_si(vektoro vkt, signed long int x, signed long int y, signed long int z);
inline void vektoron_realan_algxustigu_si(vektoro_reala vkt, signed long int x1, unsigned long int x2, signed long int y1, unsigned long int y2, signed long int z1, unsigned long int z2);

inline void vektoro_kopiu(vektoro vkt1, vektoro vkt2);
inline void vektoro_reala_kopiu(vektoro_reala vkt1, vektoro_reala vkt2);

inline void vektoron_deprenu(vektoro finvkt, vektoro vkt1, vektoro vkt2);
inline void vektoron_realan_deprenu(vektoro_reala finvkt, vektoro_reala vkt1, vektoro_reala vkt2);
inline void vektoron_aldonu(vektoro finvkt, vektoro vkt1, vektoro vkt2);
inline void vektoron_realan_aldonu(vektoro_reala finvkt, vektoro_reala vkt1, vektoro_reala vkt2);

inline void tmpvariablojn_init(struct tmpvariabloj *tmpvariabloj);
inline void tmpvariablojn_klarigu(struct tmpvariabloj *tmpvariabloj);

//void kolisiojn_rekomencigu(kolisio *kolisioj);

void longo_de_vektoro_vktmp(mpz_t longo, vektoro vkt, vektoro vktmp);
void longo_de_vektoro_reala_vktmp(mpf_t longo, vektoro_reala vkt, vektoro_reala vktmp);

void vektoron_mul_per_vektoronumbroj(vektoro vkt_fina, vektoro vkt1, vektoro vkt2);
void vektoron_realan_mul_per_vektoronumbroj(vektoro_reala vkt_fina, vektoro_reala vkt1, vektoro_reala vkt2);
inline void vektoron_realan_mul_per_numbro(vektoro_reala vkt_fina, vektoro_reala vkt, mpq_t numbro);

void vektoron_dot(mpz_t dot, vektoro vkt1, vektoro vkt2, vektoro vktmp);
void vektoron_realan_dot(mpq_t dot, vektoro_reala vkt1, vektoro_reala vkt2);

void vektoron_dividu_per_mpz(vektoro vkt, mpz_t numbro);
void vektoron_realan_dividu_per_mpq(vektoro_reala vkt, mpq_t *numbro);

void vektoro_de_unuo(vektoro vkt, vektoro vkt_du);
void vektoro_reala_de_unuo(vektoro_reala vkt, vektoro_reala vkt_du);

void montru_al_loko(vektoro loko, duonrekto dr, mpz_t t);
void montru_al_loko_qq(vektoro_reala loko, duonrektoqq dr, mpq_t t);

void duonrektoqq_tusxas_sferon(mpf_t respondo, vektoro_reala centrodesfero, mpq_t radiuso, duonrektoqq dr, vektoro_reala lmc, struct tmpvariabloj *tmpvariabloj, struct konstantoj *konstantoj);

void duonrekto_tusxas_sferon(mpz_t respondo, struct sfero *sfero, duonrekto dr, struct tmpvariabloj *tv, struct konstantoj *konst);

//struct sfero {
//  vektoro centro;
//  mpz_t radiuso;
//  materialo *materialo;
//};

void duonrekto_tusxas_sferon(mpz_t respondo, struct sfero *sfero, duonrekto dr, struct tmpvariabloj *tv, struct konstantoj *konst){
  //lmc = dr[0]-centrodesfero
  vektoron_deprenu(tv->vktmp, dr[0], sfero->centro);
  //a = dot(dr[1], dr[1]);
  vektoron_dot(tv->vktmp2[0], dr[1], dr[1]);
  //b = 2*dot(lmc, dr[1]);
  vektoron_dot(tv->vktmp2[1], tv->vktmp, dr[1]);
  mpz_add(tv->vktmp2[1], tv->vktmp2[1], tv->vktmp2[1]);
  //c = dot(lmc, lmc) - radiuso*radiuso;
  mpz_mul(tv->mpztmp, radiuso, radiuso);
  vektoron_dot(tv->vktmp2[2], tv->vktmp, tv->vktmp);
  mpz_sub(tv->vktmp2[2], tv->vktmp[2], tv->mpztmp);
  //diskriminento = b*b - 4*a*c;
  
  //diskriminento estas mpqtmp

  //if (diskriminento<0){return -1;}

  //else {return (-b - sqrt(discriminant) ) / (2.0*a);}
}

void duonrektoqq_tusxas_sferon(mpf_t respondo, vektoro_reala centrodesfero, mpq_t radiuso, duonrektoqq dr, vektoro_reala lmc, struct tmpvariabloj *tv, struct konstantoj *konstantoj){
  //lmc = dr[0]-centrodesfero
  vektoron_realan_deprenu(lmc, dr[0], centrodesfero);
  //a = dot(dr[1], dr[1]);
  vektoron_realan_dot(tv->vkrtmp[0], dr[1], dr[1]);
  //b = 2*dot(lmc, dr[1]);
  vektoron_realan_dot(tv->vkrtmp[1], lmc, dr[1]);
  mpq_mul(tv->vkrtmp[1], tv->vkrtmp[1], konstantoj->du);
  //c = dot(lmc, lmc) - radiuso*radiuso;
  vektoron_realan_dot(tv->vkrtmp[2], lmc, lmc);
  mpq_mul(tv->mpqtmp, radiuso, radiuso);
  mpq_sub(tv->vkrtmp[2], tv->vkrtmp[2], tv->mpqtmp);
  //diskriminento = b*b - 4*a*c;
  mpq_mul(tv->mpqtmp2, tv->vkrtmp[1], tv->vkrtmp[1]);
  mpq_mul(tv->mpqtmp, tv->vkrtmp[0], konstantoj->kvar);
  mpq_mul(tv->mpqtmp, tv->mpqtmp, tv->vkrtmp[2]);
  mpq_sub(tv->mpqtmp, tv->mpqtmp2, tv->mpqtmp);
  //diskriminento estas mpqtmp

  //malnova
  ////return (diskriminento>0);
  //return (mpq_cmp_ui(mpqtmp, 0, 1)>0);

  //nova
  //if (diskriminento<0){return -1;}
  if (mpq_cmp_si(tv->mpqtmp, 0, 1)<0){
    mpf_set_d(respondo, -1);
  }
  //else {return (-b - sqrt(discriminant) ) / (2.0*a);}
  else {
    mpf_set_q(respondo, tv->mpqtmp);
    mpf_sqrt(respondo, respondo);
    mpf_set_q(tv->mpftmp, tv->vkrtmp[1]);
    mpf_neg(tv->mpftmp, tv->mpftmp);
    mpf_sub(respondo, tv->mpftmp, respondo);

    mpf_set_q(tv->mpftmp, tv->vkrtmp[0]);
    mpf_mul_ui(tv->mpftmp, tv->mpftmp, 2);
    mpf_div(respondo, respondo, tv->mpftmp);
  }
}

void montru_al_loko_qq(vektoro_reala loko, duonrektoqq dr, mpq_t t){
  //loko = dr[0] + t*dr[1];
  vektoron_realan_mul_per_numbro(loko, dr[1], t);
  vektoron_realan_aldonu(loko, loko, dr[0]);
}

void vektoron_realan_deprenu(vektoro_reala finvkt, vektoro_reala vkt1, vektoro_reala vkt2){
  mpq_sub(finvkt[0], vkt1[0], vkt2[0]);
  mpq_sub(finvkt[1], vkt1[1], vkt2[1]);
  mpq_sub(finvkt[2], vkt1[2], vkt2[2]);
}
void vektoron_deprenu(vektoro finvkt, vektoro vkt1, vektoro vkt2){
  mpz_sub(finvkt[0], vkt1[0], vkt2[0]);
  mpz_sub(finvkt[1], vkt1[1], vkt2[1]);
  mpz_sub(finvkt[2], vkt1[2], vkt2[2]);
}
void vektoron_realan_aldonu(vektoro_reala finvkt, vektoro_reala vkt1, vektoro_reala vkt2){
  mpq_add(finvkt[0], vkt1[0], vkt2[0]);
  mpq_add(finvkt[1], vkt1[1], vkt2[1]);
  mpq_add(finvkt[2], vkt1[2], vkt2[2]);
}
void vektoron_aldonu(vektoro finvkt, vektoro vkt1, vektoro vkt2){
  mpz_add(finvkt[0], vkt1[0], vkt2[0]);
  mpz_add(finvkt[1], vkt1[1], vkt2[1]);
  mpz_add(finvkt[2], vkt1[2], vkt2[2]);
}

void longo_de_vektoro_reala_vktmp(mpf_t longo, vektoro_reala vkt, vektoro_reala vktmp){
  vektoron_realan_mul_per_vektoronumbroj(vktmp, vkt, vkt);
  //aldonu ĉiujn numbrojn en vktmp al vktmp[0]
  mpq_add(vktmp[0], vktmp[0], vktmp[1]);
  mpq_add(vktmp[0], vktmp[0], vktmp[2]);

  mpf_set_q(longo, vktmp[0]);
  //sqrt longon
  mpf_sqrt(longo, longo);
}

void vektoron_realan_dot(mpq_t dot, vektoro_reala vkt1, vektoro_reala vkt2){
  vektoro_reala vktmp;
  vektoron_realan_init(vktmp);
  vektoron_realan_mul_per_vektoronumbroj(vktmp, vkt1, vkt2);
  mpq_add(dot, vktmp[0], vktmp[1]);
  mpq_add(dot, dot, vktmp[2]);
  vektoron_realan_klarigu(vktmp);
}

void vektoron_dot(mpz_t dot, vektoro vkt1, vektoro vkt2, vektoro vktmp){
  vektoron_mul_per_vektoronumbroj(vktmp, vkt1, vkt2);
  mpz_add(dot, vktmp[0], vktmp[1]);
  mpz_add(dot, dot, vktmp[2]);
}

void vektoron_realan_mul_per_vektoronumbroj(vektoro_reala vkt_fina, vektoro_reala vkt1, vektoro_reala vkt2){
  mpq_mul(vkt_fina[0], vkt1[0], vkt2[0]);
  mpq_mul(vkt_fina[1], vkt1[1], vkt2[1]);
  mpq_mul(vkt_fina[2], vkt1[2], vkt2[2]);
}

void vektoron_realan_mul_per_numbro(vektoro_reala vkt_fina, vektoro_reala vkt, mpq_t numbro){
  mpq_mul(vkt_fina[0], vkt[0], numbro);
  mpq_mul(vkt_fina[1], vkt[1], numbro);
  mpq_mul(vkt_fina[2], vkt[2], numbro);
}

void duonrektonqq_init(duonrektoqq dr){
  vektoron_realan_init(dr[0]);
  vektoron_realan_init(dr[1]);
}
void vektoron_realan_init(vektoro_reala vkt){
  mpq_init(vkt[0]);
  mpq_init(vkt[1]);
  mpq_init(vkt[2]);
}
void vektoron_init(vektoro vkt){
  mpz_init(vkt[0]);
  mpz_init(vkt[1]);
  mpz_init(vkt[2]);
}

//void kolisiojn_rekomencigu(kolisio *kolisioj){
//  for (j=0; j<(kolisioj->numbro_da_kolisiovektoroj); j++){
//    vektoron_klarigu(kolisioj->kolisio[j].vektoro);
//  }
//  free(kolisioj->kolisio);
//  kolisioj->kolisio = malloc(sizeof(struct kolisiovektoro));
//  vektoron_init(kolisioj->kolisio[0].vektoro);
//}

//struct kolisiovektoro {
//  vektoro_reala vektoro;
//  materialo *materialo;
//} kolisiovektoro;
//typedef struct kolisioj {
//  unsigned char estas;
//  vektoro kolisiejo;
//  size_t numbro_da_kolisiovektoroj;
//  kolisiovektoro *kolisio;
//} kolisio;

void duonrektonqq_klarigu(duonrektoqq dr){
  vektoron_realan_klarigu(dr[0]);
  vektoron_realan_klarigu(dr[1]);
}
void vektoron_realan_klarigu(vektoro_reala vkt){
  mpq_clear(vkt[0]);
  mpq_clear(vkt[1]);
  mpq_clear(vkt[2]);
}
void vektoron_klarigu(vektoro vkt){
  mpz_clear(vkt[0]);
  mpz_clear(vkt[1]);
  mpz_clear(vkt[2]);
}

void vektoro_reala_kopiu(vektoro_reala vkt1, vektoro_reala vkt2){
  mpq_set(vkt1[0], vkt2[0]);
  mpq_set(vkt1[1], vkt2[1]);
  mpq_set(vkt1[2], vkt2[2]);
}
void vektoro_kopiu(vektoro vkt1, vektoro vkt2){
  mpz_set(vkt1[0], vkt2[0]);
  mpz_set(vkt1[1], vkt2[1]);
  mpz_set(vkt1[2], vkt2[2]);
}

void vektoron_realan_algxustigu(vektoro_reala vkt, mpq_t x, mpq_t y, mpq_t z){
  mpq_set(vkt[0], x);
  mpq_set(vkt[1], y);
  mpq_set(vkt[2], z);
}
void vektoron_algxustigu(vektoro vkt, mpz_t x, mpz_t y, mpz_t z){
  mpz_set(vkt[0], x);
  mpz_set(vkt[1], y);
  mpz_set(vkt[2], z);
}

void vektoron_realan_algxustigu_si(vektoro_reala vkt, signed long int x1, unsigned long int x2, signed long int y1, unsigned long int y2, signed long int z1, unsigned long int z2){
  mpq_set_si(vkt[0], x1, x2);
  mpq_set_si(vkt[1], y1, y2);
  mpq_set_si(vkt[2], z1, z2);
}
void vektoron_algxustigu_si(vektoro vkt, signed long int x, signed long int y, signed long int z){
  mpz_set_si(vkt[0], x);
  mpz_set_si(vkt[1], y);
  mpz_set_si(vkt[2], z);
}

void tmpvariablojn_init(struct tmpvariabloj *tmpvariabloj){
  mpq_init(tmpvariabloj->mpqtmp);
  mpq_init(tmpvariabloj->mpqtmp2);
  mpz_init(tmpvariabloj->mpztmp);
  mpz_init(tmpvariabloj->mpztmp2);
  mpf_init(tmpvariabloj->mpftmp);
  vektoron_realan_init(tmpvariabloj->vkrtmp);
  vektoron_realan_init(tmpvariabloj->vkrtmp2);
}
void tmpvariablojn_klarigu(struct tmpvariabloj *tmpvariabloj){
  mpq_clear(tmpvariabloj->mpqtmp);
  mpq_clear(tmpvariabloj->mpqtmp2);
  mpz_clear(tmpvariabloj->mpztmp);
  mpz_clear(tmpvariabloj->mpztmp2);
  mpf_clear(tmpvariabloj->mpftmp);
  vektoron_realan_klarigu(tmpvariabloj->vkrtmp);
  vektoron_realan_klarigu(tmpvariabloj->vkrtmp2);
}

#endif
