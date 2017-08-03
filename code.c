#include <stdio.h>
#include <string.h>

char *show(const Galton *G) {
    int len = 15 + 10*5 + 6 + 9 + 13*3;
    char *out = (char *)malloc(len + 30*G->nsph + 1);
    double L[2][2];

    d2mv(L[0], G->frame, G->L[0]);
    d2mv(L[1], G->frame, G->L[1]);

    snprintf(out, len+1, "Galton object:\n"
                   " gamma = %6.2f\n"
                   " sigma = %8.3e\n"
                   "     g = %6.2f %6.2f\n"
                   "     L = %6.2f %6.2f\n"
                   "         %6.2f %6.2f\n",
                G->gamma, G->sigma,
                G->g*G->frame[0][0], G->g*G->frame[1][0],
                L[0][0], L[0][1],
                L[1][0], L[1][1]);

    char *buf = out + len;
    const double *v = G->sph_crds;
    for(int i=0; i<G->nsph; i++, v += 3) {
        d2mv(L[0], G->frame, v);
        snprintf(buf, 31, "   sph = %6.2f %6.2f %6.2f\n",
                            L[0][0], L[0][1], G->sph_crds[2]);
        buf += 30;
    }

    return out;
}

size_t size(const Galton *G) {
    size_t len = sizeof(double)*GALTON_DOUBLES;
    len += size_uint32(G->nsph);
    len += sizeof(double)*3*G->nsph;
    return len;
}

void serialize(SWriter *s, const Galton *G) {
    write_uint32(s, G->nsph);
    const double *v = (const double *)G;
    for(int i=0; i<GALTON_DOUBLES; i++, v++) {
        write_fixed64(s, *v);
    }
    v = G->sph_crds;
    for(int i=0; i<3*G->nsph; i++, v++) {
        write_fixed64(s, *v);
    }
}

void parse(sil_State *S, const uint8_t *buf, size_t len) {
    uint32_t nsph;
    unsigned k = read_uint32(&nsph, buf, len);
    buf += k; len -= k;

    if(len != sizeof(double)*(GALTON_DOUBLES + 3*nsph)) {
        sil_err(S, "Invalid Galton length!");
        return;
    }
    Galton *G = (Galton *)malloc(sizeof(Galton) + sizeof(double)*3*nsph);
    double *v = (double *)G;
    for(int i=0; i<GALTON_DOUBLES; i++, v++) {
        k = read_fixed64(v, buf, len);
        buf += k; len -= k;
    }
    G->co_linear = G->L[0][1] == 0.0;
    G->nsph = nsph;

    v = G->sph_crds;
    for(int i=0; i<3*nsph; i++) {
        k = read_fixed64(v, buf, len);
        buf += k; len -= k;
    }
    sil_newuserdata(S, galton_hash, G);
}

void copy(sil_State *S) {
    size_t len;
    Galton *G = (Galton *)sil_getST(S, &len);
    if(G == NULL) {
        sil_err(S, "Can't copy - no object present.");
        return;
    }
    len = sizeof(Galton) + sizeof(double)*3*G->nsph;
    Galton *Gp = (Galton *)malloc(len);
    memcpy(Gp, G, len);
    sil_setST(S, Gp, len);
}

void handler(Galton *G) {
    free(G);
}
