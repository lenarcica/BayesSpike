
/* Flags */
struct sxpinfo_struct2 {
    SEXPTYPE type      :  5;/* ==> (FUNSXP == 99) %% 2^5 == 3 == CLOSXP
			     * -> warning: `type' is narrower than values
			     *              of its type
			     * when SEXPTYPE was an enum */
    unsigned int obj   :  1;
    unsigned int named :  2;
    unsigned int gp    : 16;
    unsigned int mark  :  1;
    unsigned int debug :  1;
    unsigned int trace :  1;  /* functions and memory tracing */
    unsigned int spare :  1;  /* currently unused */
    unsigned int gcgen :  1;  /* old generation number */
    unsigned int gccls :  3;  /* node class */
}; /*		    Tot: 32 */

struct vecsxp_struct2 {
    R_len_t	length;
    R_len_t	truelength;
};

struct primsxp_struct2 {
    int offset;
};

struct symsxp_struct2 {
    struct SEXPREC *pname;
    struct SEXPREC *value;
    struct SEXPREC *internal;
};

struct listsxp_struct2 {
    struct SEXPREC *carval;
    struct SEXPREC *cdrval;
    struct SEXPREC *tagval;
};

struct envsxp_struct2 {
    struct SEXPREC *frame;
    struct SEXPREC *enclos;
    struct SEXPREC *hashtab;
};

struct closxp_struct2 {
    struct SEXPREC *formals;
    struct SEXPREC *body;
    struct SEXPREC *env;
};

struct promsxp_struct2 {
    struct SEXPREC *value;
    struct SEXPREC *expr;
    struct SEXPREC *env;
};

typedef struct SEXPREC2 {
    ;
    union {
	struct primsxp_struct2 primsxp;
	struct symsxp_struct2 symsxp;
	struct listsxp_struct2 listsxp;
	struct envsxp_struct2 envsxp;
	struct closxp_struct2 closxp;
	struct promsxp_struct2 promsxp;
    } u;
} SEXPREC2, *SEXP2;

#define SEXPREC_HEADER2 \
    struct sxpinfo_struct2 sxpinfo; \
    struct SEXPREC *attrib; \
    struct SEXPREC *gengc_next_node, *gengc_prev_node
    
#ifndef AVECSEXP2 
typedef struct VECTOR_SEXPREC2 {
    SEXPREC_HEADER2;
    struct vecsxp_struct2 vecsxp;
} VECTOR_SEXPREC2, *VECSEXP2;
#define AVECSEXP2 0
#endif
