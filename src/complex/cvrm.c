/*
    module: cvrm.c
    Purpose: miscellaneous cover manipulation
	a) verify two covers are equal, check consistency of a cover
	b) unravel a multiple-valued cover into minterms
	c) sort covers
*/

#include "cvr.h"

static pcover Fmin;
static pcube phase;

/*
 *  minimize each output function individually
 */
void so_espresso(pPLA PLA, int strategy)
{
    Fmin = new_cover(PLA->F->count);
    if (strategy == 0) {
	foreach_output_function(PLA, so_do_espresso, so_save);
    } else {
	foreach_output_function(PLA, so_do_exact, so_save);
    }
    sf_free(PLA->F);
    PLA->F = Fmin;
}


/*
 *  minimize each output function, choose function or complement based on the
 *  one with the fewer number of terms
 */
void so_both_espresso(pPLA PLA, int strategy)
{
    phase = set_save(cube.fullset);
    Fmin = new_cover(PLA->F->count);
    if (strategy == 0) {
	foreach_output_function(PLA, so_both_do_espresso, so_both_save);
    } else {
	foreach_output_function(PLA, so_both_do_exact, so_both_save);
    }
    sf_free(PLA->F);
    PLA->F = Fmin;
    PLA->phase = phase;
}


int so_do_espresso(pPLA PLA, int i)
{
    char word[32];

    /* minimize the single-output function (on-set) */
    skip_make_sparse = 1;
    (void) sprintf(word, "ESPRESSO-POS(%d)", i);
    EXEC_S(PLA->F = espresso(PLA->F, PLA->D, PLA->R), word, PLA->F);
    return 1;
}


int so_do_exact(pPLA PLA, int i)
{
    char word[32];

    /* minimize the single-output function (on-set) */
    skip_make_sparse = 1;
    (void) sprintf(word, "EXACT-POS(%d)", i);
    EXEC_S(PLA->F = minimize_exact(PLA->F, PLA->D, PLA->R, 1), word, PLA->F);
    return 1;
}


/*ARGSUSED*/
int so_save(pPLA PLA, int i)
{
    Fmin = sf_append(Fmin, PLA->F);	/* disposes of PLA->F */
    PLA->F = NULL;
    return 1;
}


int so_both_do_espresso(pPLA PLA, int i)
{
    char word[32];

    /* minimize the single-output function (on-set) */
    (void) sprintf(word, "ESPRESSO-POS(%d)", i);
    skip_make_sparse = 1;
    EXEC_S(PLA->F = espresso(PLA->F, PLA->D, PLA->R), word, PLA->F);

    /* minimize the single-output function (off-set) */
    (void) sprintf(word, "ESPRESSO-NEG(%d)", i);
    skip_make_sparse = 1;
    EXEC_S(PLA->R = espresso(PLA->R, PLA->D, PLA->F), word, PLA->R);

    return 1;
}


int so_both_do_exact(pPLA PLA, int i)
{
    char word[32];

    /* minimize the single-output function (on-set) */
    (void) sprintf(word, "EXACT-POS(%d)", i);
    skip_make_sparse = 1;
    EXEC_S(PLA->F = minimize_exact(PLA->F, PLA->D, PLA->R, 1), word, PLA->F);

    /* minimize the single-output function (off-set) */
    (void) sprintf(word, "EXACT-NEG(%d)", i);
    skip_make_sparse = 1;
    EXEC_S(PLA->R = minimize_exact(PLA->R, PLA->D, PLA->F, 1), word, PLA->R);

    return 1;
}


int so_both_save(pPLA PLA, int i)
{
    if (PLA->F->count > PLA->R->count) {
	sf_free(PLA->F);
	PLA->F = PLA->R;
	PLA->R = NULL;
	i += cube.first_part[cube.output];
	set_remove(phase, i);
    } else {
	sf_free(PLA->R);
	PLA->R = NULL;
    }
    Fmin = sf_append(Fmin, PLA->F);
    PLA->F = NULL;
    return 1;
}
